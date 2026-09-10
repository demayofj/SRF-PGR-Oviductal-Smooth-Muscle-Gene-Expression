"""
MUSCLEMOTION 'Contraction' trace analysis pipeline
==================================================
Oviduct whole-tissue video recordings, ~8.5 fps, ~10 min each.
Ryan Marquardt 2026 revised 9/10/26
NIEHS DeMayo Lab
Developed with HHS Claude Opus 4.8

WHAT THIS DOES, IN PLAIN TERMS
------------------------------
MUSCLEMOTION gives you, per recording, a text file with two columns:
    column 1 = time (ms)   <-- UNRELIABLE if framerate was rounded (MUSCLEMOTION does this by default but should be fixed now with revised macro - see below)
    column 2 = 'Contraction' in arbitrary units (a.u.)
The a.u. value is a whole-field displacement-from-reference measure. It is NOT
force or diameter, its absolute size is NOT comparable between recordings, and
these traces carry a large slow BASELINE DRIFT that must be removed before any
amplitude or peak analysis means anything.

This script turns each trace into a handful of scalar summary numbers
(frequency, interval, regularity, duty cycle, within-recording ratio of mean to median peak height), ONE ROW
PER RECORDING, saved to a CSV. Those per-recording rows are what you then
compare between genotypes -- with ANIMAL (not contraction) as the unit of n.

It also writes a QC plot per file so you can eyeball whether the detrending and
peak-calling behaved. LOOK AT THESE. They are the analogue of the MaxProjectStack
QC check -- they are how you catch a bad detrend window or mis-tuned threshold.

HOW TO RUN
----------
Interactively (stepping through), just make sure INPUT_DIR below points at your
folder, then run blocks in order, or run the whole file.

Command line checks to ensure python is ready to go:
    python3 --version   [I have Python 3.10.8 installed on my NIEHS HPC user]
    which python3   [~/miniforge3/bin/python3]
    readlink -f $(which python3)    [/ddn/gs1/home/marquardtrm/miniforge3/bin/python3.10]


From the command line / a job script:
    export MPLBACKEND=Agg          # headless: save plots to file, never display
    To run interactively, enter the interpreter and proceed below: python3
    To run full script: python3 [this-script.py] --input-dir [insert-input-directory]

Requires: numpy, pandas, scipy, matplotlib  (already present in your miniforge base)
"""

import os
import glob
import argparse
import numpy as np
import pandas as pd
from scipy.signal import find_peaks, savgol_filter
import matplotlib
# Force the non-interactive 'Agg' backend BEFORE importing pyplot. On a headless
# compute node there is no display; Agg saves figures to disk and never tries to
# open a window. Setting the MPLBACKEND env var does the same thing; doing it here
# too means the script is safe even if you forget the env var.
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ==========================================================================
# CONFIG  --  the knobs you are most likely to touch
# ==========================================================================

# Where your MUSCLEMOTION 'contraction' .txt files live. Can be overridden on
# the command line with --input-dir. When stepping through interactively, just
# edit this string to your folder.
INPUT_DIR = "."

# Glob pattern used to FIND the contraction files inside INPUT_DIR. MUSCLEMOTION
# names them 'contraction.txt' (and 'contraction1.txt', etc. on repeat runs), so
# this pattern grabs all of them while ignoring the 'speed-of-contraction' files.
# Adjust if your naming differs.
FILE_GLOB = "*contraction*.txt"

# Files whose names contain any of these substrings are EXCLUDED. We drop the
# speed files, which also contain the word 'contraction'.
EXCLUDE_SUBSTRINGS = ["speed"]

FPS = 8.5                      # TRUE recording framerate. The macro may have
                               # stored 9 (its dialog rounded 8.5), so we do NOT
                               # trust the file's time column -- we rebuild the
                               # time axis from this value. Set it correctly.

TRIM_START_S = 0               # Discard the first N seconds. The very start of a
                               # MUSCLEMOTION trace often has settling / reference-
                               # frame artifacts.

DETREND_WINDOW_S = 120.0       # Rolling-baseline window (seconds). MUST be much
                               # longer than one contraction period, or it will
                               # chew into the peaks. Your periods look ~30-70 s,
                               # so 120 s is a safe start. CHECK on the QC
                               # plot: the baseline line should trace the slow
                               # drift, NOT wiggle with each contraction.

MIN_PEAK_DISTANCE_S = 1.0      # Refractory period: two peaks cannot be counted
                               # closer than this. Stops noise wiggles near a
                               # crest from being counted as separate contractions.
                               # NOTE: 1.0 s is permissive for SLOW oviduct
                               # contractions -- it will catch small/fast events
                               # but also risks counting noise near a crest as
                               # multiple peaks. VERIFY against the QC plots (and
                               # ideally the video). Raise to 5-12 s if you see
                               # doublets on single contractions.

SMOOTH_WINDOW_S = 1.0          # Light smoothing before peak detection, to kill
                               # single-frame jitter without blurring real peaks.
                               # With MIN_PEAK_DISTANCE_S also low, the detector
                               # is sensitive to small fast fluctuations -- lean
                               # on the QC plots to confirm peaks are real.

PROMINENCE_SD_MULT = 0.1       # Peak height gate = this multiple of the trace's
                               # own robust noise level (see detect_peaks). Raise
                               # to be stricter, lower to catch smaller events.
                               # TUNE THIS AGAINST THE QC PLOT.

#  --- GENOTYPE ASSIGNMENT BY FILENAME KEYWORD --------------------------------
# Genotype is detected from a substring in the filename (case-insensitive).
# In this case, filenames contain either "srfflox" (control) or "PRcreSrf" (mutant).
# Update as appropriate.
#
# ORDER MATTERS: first match wins. "prcresrf" is checked before "srfflox" as a
# safeguard, in case a mutant name ever also contains an "srf..." fragment.
GENOTYPE_KEYWORDS = [
    ("prcresrf", "Pgr_cre_Srf_ff"),   # mutant
    ("srfflox",  "Srf_ff"),           # control
]


# ==========================================================================
# STAGE 1 -- find the files
# ==========================================================================
def discover_files(input_dir):
    """Return a sorted list of contraction .txt paths in input_dir,
    excluding speed files. This replaces hand-listing every filename."""
    pattern = os.path.join(input_dir, FILE_GLOB)
    hits = sorted(glob.glob(pattern))
    # Drop anything matching an exclude substring (e.g. 'speed-of-contraction')
    kept = [f for f in hits
            if not any(bad in os.path.basename(f).lower()
                       for bad in EXCLUDE_SUBSTRINGS)]
    return kept


# ==========================================================================
# STAGE 2 -- load one trace
# ==========================================================================
def load_trace(path):
    """Read the two-column file and REBUILD the time axis from FPS.

    Why rebuild instead of using column 1? Because MUSCLEMOTION's dialog only
    accepts an integer framerate by default, so 8.5 fps may have been stored as 9, which
    would stretches the file's time column by ~5.9%. The CONTRACTION VALUES and the
    frame ORDER are unaffected, so we keep column 2 and regenerate time as
    (frame index) * (1000 / FPS). The MUSUCLEMOTION macro can be revised to fix this
    (which I did in this case, but it doesn't hurt to rebuild it anyway here to be safe"""
    df = pd.read_csv(path, sep=r"\s+", header=None,
                     names=["time_ms_file", "au"], engine="python")
    df = df.dropna().reset_index(drop=True)

    dt = 1000.0 / FPS                         # true sampling interval in ms
    df["time_ms"] = np.arange(len(df)) * dt   # clean, correct time axis

    # If the file's own interval disagrees with the true one, say so once.
    # (Confirms the rounding happened and that we were right to rebuild.)
    if len(df) > 1:
        file_dt = np.median(np.diff(df["time_ms_file"].values))
        if abs(file_dt - dt) / dt > 0.005:
            print(f"    [note] file interval {file_dt:.2f} ms vs true {dt:.2f} ms "
                  f"-> rebuilt time axis from FPS={FPS}")
    return df, dt


# ==========================================================================
# STAGE 3 -- clean: trim the start, remove zero-dropout artifacts
# ==========================================================================
def clean_trace(df):
    """Trim the first TRIM_START_S seconds (if needed) and interpolate over any 0 values,
    which MUSCLEMOTION occasionally emits as dropouts (not real contractions)."""
    df = df[df["time_ms"] >= TRIM_START_S * 1000.0].copy()
    df.loc[df["au"] == 0, "au"] = np.nan          # mark zeros as missing
    df["au"] = df["au"].interpolate(limit_direction="both")  # fill them in
    return df.reset_index(drop=True)


# ==========================================================================
# STAGE 4 -- DETREND
# ==========================================================================
def detrend(df, dt_ms):
    """Estimate the slow baseline drift with a rolling MEDIAN and subtract it,
    leaving the oscillations centered near zero.
    
    Because MUSLCEMOTION uses a fixed early reference frame, a monotonic rise
    in the raw output is normal due to the tissue progressively moving/displacing 
    away from its frame-1 position — slow drift, swelling, settling, or focus change,
    not a real rise in contraction amplitude. Detrending levels this drift.
    
    As long as the rolling median window spans well over one
    contraction period, the median tracks the slow drift and ignores the
    (comparatively brief) contraction humps. Subtracting it flattens the
    baseline so that peak heights and thresholds finally mean the same thing
    across the whole recording. Check this on the QC outputs and adjust as needed"""
    n_win = int(round(DETREND_WINDOW_S * 1000.0 / dt_ms))  # window in samples
    if n_win % 2 == 0:                                      # make it odd so it
        n_win += 1                                          # can be centered
    au = df["au"].values.astype(float)

    baseline = (pd.Series(au)
                .rolling(n_win, center=True, min_periods=1)
                .median()
                .values)
    df["baseline"] = baseline        # keep for the QC plot
    df["detr"] = au - baseline        # detrended signal: oscillations about 0
    return df


# ==========================================================================
# STAGE 5 -- smooth lightly, then detect peaks (= contractions)
# ==========================================================================
def detect_peaks(df, dt_ms):
    """Lightly smooth the detrended trace, then find peaks using two guards:
      * a minimum inter-peak DISTANCE (refractory period), and
      * a PROMINENCE threshold set from the trace's own noise level.
    Prominence (how far a peak stands out from the surrounding trace) is a
    better gate than raw height on a detrended signal."""
    x = df["detr"].values

    # Savitzky-Golay smoothing: fits a small polynomial in a sliding window,
    # which suppresses jitter while preserving peak shape better than a moving
    # average. Window must be odd and >= polyorder+2.
    w = int(round(SMOOTH_WINDOW_S * 1000.0 / dt_ms))
    if w % 2 == 0:
        w += 1
    w = max(w, 5)
    xs = savgol_filter(x, window_length=w, polyorder=2)
    df["smooth"] = xs

    # Minimum spacing between peaks, in samples.
    min_dist = int(round(MIN_PEAK_DISTANCE_S * 1000.0 / dt_ms))
    min_dist = max(min_dist, 1)     # find_peaks requires distance >= 1

    # Robust noise estimate: the median absolute deviation (MAD) scaled to be
    # comparable to a standard deviation. Robust = not thrown off by the peaks
    # themselves, unlike a plain std. The peak must stand out by at least
    # PROMINENCE_SD_MULT times this noise level.
    noise = np.median(np.abs(xs - np.median(xs))) / 0.6745
    prominence = PROMINENCE_SD_MULT * noise

    peaks, props = find_peaks(xs, distance=min_dist, prominence=prominence)
    return peaks, props, noise


# ==========================================================================
# STAGE 6 -- reduce one recording to scalar summary metrics
# ==========================================================================
def per_minute_counts(df, peaks, n_minutes=10):
    """Count detected contractions in each 1-minute bin of the recording.

    Bins are defined on the REBUILT (FPS-corrected) time axis in absolute
    recording time: minute 1 = [0,60) s, minute 2 = [60,120) s, etc. This makes
    the minute columns line up across files regardless of TRIM_START_S.

    Returns a dict {'contractions_min_1': n1, ..., 'contractions_min_10': n10}.
    A minute that lies entirely BEYOND the end of the recording is reported as
    NaN (not 0), so a short recording doesn't masquerade as zero activity in
    minutes it never captured. A minute that is only PARTIALLY covered still
    reports its count (it will just be over less than a full 60 s)."""
    pk_times_s = df["time_ms"].values[peaks] / 1000.0    # peak times in seconds
    rec_end_s = df["time_ms"].iloc[-1] / 1000.0          # last sample's time

    out = {}
    for minute in range(1, n_minutes + 1):
        lo = (minute - 1) * 60.0
        hi = minute * 60.0
        col = f"contractions_min_{minute}"
        if lo >= rec_end_s:
            # This whole minute is past the end of the recording -> not observed.
            out[col] = np.nan
        else:
            # Count peaks whose time falls in [lo, hi).
            out[col] = int(np.sum((pk_times_s >= lo) & (pk_times_s < hi)))
    return out


def summarize(df, peaks, dt_ms):
    """Collapse the whole trace into a few numbers. Ordered here roughly by how
    much they can be trusted at 8.5 fps (most trustworthy first)."""
    dur_s = (df["time_ms"].iloc[-1] - df["time_ms"].iloc[0]) / 1000.0
    n = len(peaks)

    # --- FREQUENCY: primary endpoint. Robust to drift and to arbitrary units. ---
    freq_per_min = n / (dur_s / 60.0)

    # --- INTERVAL and its regularity (CV) ---
    if n >= 2:
        pk_times_s = df["time_ms"].values[peaks] / 1000.0
        intervals = np.diff(pk_times_s)                 # inter-contraction gaps
        mean_ici = np.mean(intervals)
        cv_ici = np.std(intervals, ddof=1) / mean_ici    # CV: lower = more regular
    else:
        mean_ici, cv_ici = np.nan, np.nan

    # --- PEAK-HEIGHT MEAN/MEDIAN RATIO (a within-recording SHAPE descriptor,
    #     NOT an amplitude or strength measure). Because peak heights are divided
    #     by this recording's own median peak, the arbitrary-unit scaling cancels
    #     and the result reduces to mean(peak_heights)/median(peak_heights). That
    #     is a skewness indicator of the peak-height distribution (near 1.0 =
    #     symmetric; >1.0 = a few unusually large peaks). It carries NO absolute
    #     or cross-animal amplitude information -- MUSCLEMOTION a.u. cannot give
    #     that -- so do not report it as contraction strength. Raw median peak
    #     height in a.u. is kept only as a QC reference (not cross-comparable). ---
    if n >= 1:
        peak_heights = df["smooth"].values[peaks]        # heights on DETRENDED trace
        med_h = np.median(peak_heights)
        peak_mean_to_median = np.mean(peak_heights / med_h) if med_h > 0 else np.nan
        raw_amp_median = med_h
    else:
        peak_mean_to_median, raw_amp_median = np.nan, np.nan

    # --- DUTY CYCLE: fraction of time the signal sits above a low threshold. ---
    med = np.median(df["smooth"])
    mad = np.median(np.abs(df["smooth"] - med)) / 0.6745
    thr = med + mad
    duty_cycle = float(np.mean(df["smooth"].values > thr))

    result = {
        "duration_s": round(dur_s, 1),
        "n_contractions": n,
        "freq_per_min": round(freq_per_min, 3),
        "mean_interval_s": round(mean_ici, 2) if not np.isnan(mean_ici) else np.nan,
        "cv_interval": round(cv_ici, 3) if not np.isnan(cv_ici) else np.nan,
        "peak_mean_to_median": round(peak_mean_to_median, 3) if not np.isnan(peak_mean_to_median) else np.nan,
        "duty_cycle": round(duty_cycle, 3),
        "raw_amp_median_au": round(raw_amp_median, 1) if not np.isnan(raw_amp_median) else np.nan,
    }
    # Add the 10 per-minute contraction-count columns.
    result.update(per_minute_counts(df, peaks, n_minutes=10))
    return result


# ==========================================================================
# STAGE 7 -- QC plot  (LOOK AT THESE for every file)
# ==========================================================================
def qc_plot(df, peaks, animal, genotype, out_dir):
    """Two-panel figure per recording:
      top    = raw trace + the rolling baseline we subtracted
      bottom = detrended + smoothed trace with detected peaks marked
    If the baseline wiggles with each contraction, DETREND_WINDOW_S is too short.
    If peaks are missed or noise is marked, tune PROMINENCE_SD_MULT / distances."""
    fig, ax = plt.subplots(2, 1, figsize=(12, 6), sharex=True)
    t = df["time_ms"] / 1000.0

    ax[0].plot(t, df["au"], lw=0.6, label="raw")
    ax[0].plot(t, df["baseline"], lw=1.2, label="rolling baseline")
    ax[0].set_ylabel("Contraction (a.u.)")
    ax[0].legend(loc="upper left", fontsize=8)
    ax[0].set_title(f"{animal}  ({genotype})  -- raw + baseline")

    ax[1].plot(t, df["smooth"], lw=0.7, label="detrended + smoothed")
    ax[1].plot(t.values[peaks], df["smooth"].values[peaks], "rv",
               ms=6, label=f"peaks (n={len(peaks)})")
    ax[1].axhline(0, color="k", lw=0.4)
    ax[1].set_ylabel("Detrended (a.u.)")
    ax[1].set_xlabel("Time (s)")
    ax[1].legend(loc="upper left", fontsize=8)

    fig.tight_layout()
    # Sanitize the identifier for use in a filename (it's the full stem, which
    # may contain characters awkward for filenames -- keep it simple).
    safe = "".join(c if c.isalnum() or c in "-_." else "_" for c in animal)
    out_path = os.path.join(out_dir, f"QC_{safe}.png")
    fig.savefig(out_path, dpi=110)
    plt.close(fig)          # free memory -- important when looping many files
    return out_path


# ==========================================================================
# METADATA -- derive genotype from the filename
# ==========================================================================
def parse_metadata(fname):
    """Return {'genotype':..., 'animal':...} for a filename.
    Genotype comes from a case-insensitive keyword match (first match wins;
    see GENOTYPE_KEYWORDS). The full filename stem is used as the identifier --
    we do not scrape an animal number. Files with no matching keyword are
    labeled UNKNOWN (with a warning) so they still get analyzed."""
    low = fname.lower()

    genotype = "UNKNOWN"
    for keyword, label in GENOTYPE_KEYWORDS:
        if keyword in low:                # keywords in GENOTYPE_KEYWORDS are lowercase
            genotype = label
            break

    if genotype == "UNKNOWN":
        print(f"    [warn] no genotype keyword in {fname} -> UNKNOWN")

    animal = os.path.splitext(fname)[0]   # filename stem as identifier
    return {"genotype": genotype, "animal": animal}


# ==========================================================================
# STAGE 8 -- process one file end to end
# ==========================================================================
def process_file(path, out_dir):
    """Run stages 2-7 on a single file and return its summary-metrics dict."""
    fname = os.path.basename(path)

    # Derive genotype + identifier from the filename (see parse_metadata).
    meta = parse_metadata(fname)

    df, dt = load_trace(path)
    df = clean_trace(df)
    df = detrend(df, dt)
    peaks, props, noise = detect_peaks(df, dt)
    metrics = summarize(df, peaks, dt)

    metrics.update({
        "animal": meta["animal"],
        "genotype": meta["genotype"],
        "file": fname,
        "dt_ms": round(dt, 2),
        "noise_au": round(noise, 2),
    })
    qc_path = qc_plot(df, peaks, meta["animal"], meta["genotype"], out_dir)

    print(f"    {meta['animal']} ({meta['genotype']}): "
          f"n={metrics['n_contractions']}, freq={metrics['freq_per_min']}/min, "
          f"CV={metrics['cv_interval']}  ->  {os.path.basename(qc_path)}")
    return metrics


# ==========================================================================
# STAGE 9 -- driver: discover, process all, write the summary table
# ==========================================================================
def run(input_dir):
    files = discover_files(input_dir)
    if not files:
        print(f"No files matching '{FILE_GLOB}' in {input_dir} "
              f"(after excluding {EXCLUDE_SUBSTRINGS}).")
        return None

    print(f"Found {len(files)} file(s) in {input_dir}:")
    for f in files:
        print(f"    {os.path.basename(f)}")

    # Put outputs (QC plots + summary CSV) alongside the inputs.
    out_dir = input_dir
    print("\nProcessing:")
    rows = [process_file(f, out_dir) for f in files]

    summary = pd.DataFrame(rows)
    # Order columns sensibly (identity first, then metrics, then per-minute
    # counts, then diagnostics).
    per_min_cols = [f"contractions_min_{m}" for m in range(1, 11)]
    col_order = (["file", "animal", "genotype", "duration_s", "n_contractions",
                  "freq_per_min", "mean_interval_s", "cv_interval",
                  "peak_mean_to_median", "duty_cycle", "raw_amp_median_au"]
                 + per_min_cols
                 + ["dt_ms", "noise_au"])
    summary = summary[[c for c in col_order if c in summary.columns]]

    out_csv = os.path.join(out_dir, "per_recording_summary.csv")
    summary.to_csv(out_csv, index=False)
    print(f"\nWrote {out_csv}")
    print(summary.to_string(index=False))

    # ---- PRISM-READY PER-MINUTE TABLES -----------------------------------
    # For the per-minute line graph in Prism you want rows = minute (1-10) and
    # columns = animals, with a separate block per genotype. We write one CSV
    # per genotype in that orientation so each can be pasted into its Prism group.
    for geno, block in summary.groupby("genotype"):
        # take the per-minute columns, index by animal, then transpose (.T) so
        # minutes become rows and animals become columns
        t = block.set_index("animal")[per_min_cols].T
        t.index = range(1, 11)           # relabel rows 1..10 (the minutes)
        t.index.name = "minute"
        prism_path = os.path.join(out_dir, f"prism_per_minute_{geno}.csv")
        t.to_csv(prism_path)
        print(f"Wrote {prism_path}  (rows=minute, cols=animals)")

    return summary


# ==========================================================================
# ENTRY POINT
# ==========================================================================
# When run as a script this parses --input-dir. When stepping through
# interactively, this block is skipped if you import functions, OR you can just
# call run(INPUT_DIR) yourself after editing INPUT_DIR above.
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Analyze MUSCLEMOTION contraction traces.")
    parser.add_argument("--input-dir", default=INPUT_DIR,
                        help="Directory containing *contraction*.txt files "
                             "(default: current directory).")
    args = parser.parse_args()
    run(args.input_dir)