"""
================================================================================
MUSCLEMOTION 'Contraction' trace analysis pipeline  --  VERBOSE TEACHING VERSION
================================================================================
Oviduct whole-tissue video recordings, ~8.5 fps, ~10 min each.
Ryan Marquardt 2026 -- NIEHS DeMayo Lab
Developed with Claude Opus 4.8

------------------------------------------------------------------------------
WHO THIS VERSION IS FOR
------------------------------------------------------------------------------
This is a heavily-annotated copy of the analysis script. It does EXACTLY the
same thing as the normal version -- same math, same settings, same output files.
The only difference is that nearly every line is explained, aimed at a biologist
who understands programming concepts in general (variables, loops, functions,
if/else) but is newer to Python specifically. If you already read Python fluently,
use the shorter version; this one is for learning what each piece does and why.

------------------------------------------------------------------------------
THE BIG PICTURE (read this first)
------------------------------------------------------------------------------
MUSCLEMOTION produces, for each recorded oviduct video, a plain text file with
two columns:
    column 1 = time in milliseconds
    column 2 = "Contraction" in arbitrary units (a.u.)

Three facts about that 'Contraction' number drive everything this script does:

  1. It is a WHOLE-FIELD DISPLACEMENT measure (how much the tissue moved away
     from a chosen reference frame), NOT force and NOT tube diameter. Its
     absolute size is not comparable from one recording to the next, because it
     depends on tissue size, lighting, and framing.

  2. The trace has a large slow BASELINE DRIFT -- a gradual rise over the whole
     recording caused by the tissue slowly shifting away from the reference
     frame, not by real contractions getting bigger. We must REMOVE this drift
     ("detrend") before measuring anything.

  3. The real contractions are the OSCILLATIONS riding on top of that drift.
     After detrending, we detect them as PEAKS.

What the script outputs, per recording (one row per animal in a CSV):
    - how many contractions, and frequency per minute   (the main results)
    - the interval between contractions and its regularity
    - a within-recording ratio of mean to median peak height and a "duty cycle"
    - the count of contractions in each of the 10 minutes
Plus one QC (quality-control) PLOT per recording so you can SEE whether the
detrending and peak-finding did something sensible. ALWAYS look at those plots.

------------------------------------------------------------------------------
HOW A PYTHON SCRIPT LIKE THIS IS ORGANIZED (orientation for newer Python users)
------------------------------------------------------------------------------
  * "import" lines at the top load extra toolboxes (libraries) written by others.
  * ALL-CAPS names near the top are CONFIG values -- the settings you change.
  * "def name(...):" blocks define FUNCTIONS -- reusable mini-programs. Defining
    a function does NOT run it; it just teaches Python how. It runs only when
    "called" later, e.g. detrend(df, dt).
  * Indentation (the leading spaces) is how Python groups code. Everything
    indented under a "def" or "for" or "if" belongs to that block. This is not
    cosmetic in Python -- the spacing IS the structure.
  * The very bottom "if __name__ == '__main__':" block is what actually runs
    when you type  python3 thisfile.py.

------------------------------------------------------------------------------
HOW TO RUN
------------------------------------------------------------------------------
Check Python is available (type these in the terminal, not in Python):
    python3 --version            # should show 3.something
    which python3                # shows which python3 you're using

Run the whole thing on a folder of files:
    export MPLBACKEND=Agg         # tell plotting to save to file (no screen)
    python3 musclemotion-verbose.py --input-dir /path/to/your/txt/files

Or step through it interactively (see notes at the very bottom).

Requires the libraries: numpy, pandas, scipy, matplotlib
(already installed in your miniforge Python).
================================================================================
"""

# ==========================================================================
# IMPORTS -- load the toolboxes we need
# ==========================================================================
# Each "import" makes another library's tools available. We give some of them
# short nicknames with "as" so we don't have to type the full name every time.

import os        # "operating system" tools: build/join file paths, split names
import glob      # find files by a wildcard pattern (e.g. all *.txt files)
import argparse  # lets the script accept options typed on the command line
import numpy as np   # numerical arrays + math. Nickname: np
import pandas as pd  # tables ("DataFrames"), like a spreadsheet in code. Nick: pd

# From the SciPy library we pull just two specific functions (rather than all of
# it). "from LIBRARY import THING" grabs one named tool.
from scipy.signal import find_peaks, savgol_filter
#   find_peaks    -> locates peaks (local maxima) in a 1-D signal
#   savgol_filter -> a gentle smoothing filter (Savitzky-Golay)

import matplotlib
# The next line MUST come before we import matplotlib's plotting part. It tells
# matplotlib to use the "Agg" drawing engine, which writes images to files and
# never tries to open a window on screen. On a remote compute server there is no
# screen, so without this the plotting could error out or hang.
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # the actual plotting tools. Nickname: plt


# ==========================================================================
# CONFIG  --  the settings you are most likely to change
# --------------------------------------------------------------------------
# These ALL-CAPS names are "global" settings used throughout the script. By
# convention, all-caps signals "this is a constant you set once." Changing a
# value here changes the behavior everywhere that value is used.
# ==========================================================================

# Folder containing the MUSCLEMOTION 'contraction' text files. "." means "the
# current folder." You can override this on the command line with --input-dir,
# or just edit this string when stepping through interactively.
INPUT_DIR = "."

# A "glob" is a filename wildcard. The "*" matches any characters, so
# "*contraction*.txt" means: any file whose name contains the word
# "contraction" and ends in ".txt". This finds the contraction files while
# (with the exclude list below) skipping the speed files.
FILE_GLOB = "*contraction*.txt"

# Any file whose name contains one of these words is thrown out. MUSCLEMOTION's
# "speed-of-contraction" file ALSO contains the word "contraction", so we drop
# anything containing "speed". This is a Python LIST -- an ordered collection
# written in square brackets. Here it has one item.
EXCLUDE_SUBSTRINGS = ["speed"]

# The TRUE frame rate of the recordings, in frames per second. This matters
# because MUSCLEMOTION's dialog historically rounded 8.5 up to 9, which would
# stretch the time axis by ~6%. To be safe we IGNORE the file's time column and
# rebuild time ourselves from this number (see load_trace). Set it correctly.
FPS = 8.5

# Seconds to discard from the START of each recording (settling artifacts, etc.).
# 0 means "keep everything." Note: this is a number, not text.
TRIM_START_S = 0

# Detrending window in SECONDS. We estimate the slow drift with a rolling median
# over this many seconds and subtract it. This window MUST be much longer than
# one contraction period (~30-70 s here), or it will start following the
# contractions themselves and flatten them out. 120 s is a safe starting point.
# ALWAYS confirm on the QC plot that the baseline line tracks the slow drift and
# does NOT wiggle up and down with each contraction.
DETREND_WINDOW_S = 120.0

# Minimum time between two counted peaks, in seconds ("refractory period"). Two
# peaks closer than this are treated as one. 1.0 s is PERMISSIVE for these slow
# contractions -- it will catch small/fast bumps but also risks counting noise
# near a single crest as multiple contractions. VERIFY against the QC plots (and
# ideally the video). Raise to ~5-12 s if you see doubled-up peaks on one hump.
MIN_PEAK_DISTANCE_S = 1.0

# Light smoothing window in seconds, applied before peak detection to remove
# single-frame jitter. Kept short so it doesn't blur real peaks. Combined with a
# small MIN_PEAK_DISTANCE_S, the detector becomes quite sensitive to small fast
# wiggles -- again, lean on the QC plots to confirm detected peaks are real.
SMOOTH_WINDOW_S = 1.0

# The peak HEIGHT gate, expressed as a multiple of the trace's own noise level.
# A candidate must stand out ("prominence") by at least this many noise-units to
# count as a peak. Lower = catches smaller events but also more noise.
# IMPORTANT SCIENTIFIC NOTE: this is currently 0.1, which is very permissive --
# it will count very small bumps. Tune it by LOOKING AT THE QC PLOTS, never by
# adjusting until you hit a contraction count you were hoping for. Fitting the
# threshold to an expected number is how you manufacture a false result.
PROMINENCE_SD_MULT = 0.1

# --- GENOTYPE ASSIGNMENT BY FILENAME KEYWORD --------------------------------
# We read each file's genotype from a keyword in its filename. This is a LIST OF
# TUPLES. A tuple is a small fixed group of values in parentheses; here each is
# (keyword_to_look_for, genotype_label_to_assign). The keywords are lowercase
# because we compare against a lowercased filename.
#
# ORDER MATTERS: the first keyword that matches wins (see parse_metadata). We
# check the mutant keyword first as a safeguard, in case a mutant filename ever
# also contained an "srf..." fragment that could be mistaken for control.
GENOTYPE_KEYWORDS = [
    ("prcresrf", "Pgr_cre_Srf_ff"),   # if name contains "prcresrf" -> mutant
    ("srfflox",  "Srf_ff"),           # else if it contains "srfflox" -> control
]


# ==========================================================================
# STAGE 1 -- find the files
# --------------------------------------------------------------------------
# A FUNCTION named discover_files. It takes one input ("parameter") called
# input_dir and RETURNS a list of file paths. Defining it here does not run it.
# ==========================================================================
def discover_files(input_dir):
    """Return a sorted list of contraction .txt file paths inside input_dir,
    with the speed files removed.

    (Text in triple quotes right under a 'def' is a 'docstring' -- built-in
    documentation describing what the function does.)"""

    # os.path.join glues a folder and a pattern into one path string, using the
    # correct slash for this operating system. e.g. "/data" + "*contraction*.txt".
    pattern = os.path.join(input_dir, FILE_GLOB)

    # glob.glob returns every path matching the wildcard pattern. sorted(...)
    # puts them in alphabetical order so runs are reproducible.
    hits = sorted(glob.glob(pattern))

    # Now remove any file whose name contains an excluded word (like "speed").
    # This is a "list comprehension": a compact way to build a new list by
    # looping over 'hits' and keeping only items that pass a test. Read it as:
    #   "keep f, for each f in hits, IF (none of the bad words are in its name)".
    kept = [
        f
        for f in hits
        if not any(
            bad in os.path.basename(f).lower()   # is 'bad' inside the filename?
            for bad in EXCLUDE_SUBSTRINGS         # check each excluded word
        )
    ]
    #   os.path.basename(f) -> just the filename, dropping the folders
    #   .lower()            -> lowercase it, so matching ignores capitalization
    #   any(...)            -> True if ANY excluded word is present
    #   not any(...)        -> so we KEEP the file only if NO excluded word is present

    return kept   # hand the filtered list back to whoever called this function


# ==========================================================================
# STAGE 2 -- load one trace from its text file
# ==========================================================================
def load_trace(path):
    """Read one two-column file and REBUILD its time axis from FPS.

    We deliberately do NOT trust the file's own time column, because a rounded
    frame rate would have stretched it. The contraction VALUES and their ORDER
    are fine, so we keep the values and regenerate time = frame_index * (1000/FPS)
    in milliseconds."""

    # pandas reads the text file into a "DataFrame" (df) -- a table.
    #   sep=r"\s+"    -> columns are separated by whitespace (spaces/tabs)
    #   header=None   -> the file has no header row; row 1 is already data
    #   names=[...]   -> we name the two columns ourselves
    #   engine="python" -> a parser flexible enough for the \s+ separator
    df = pd.read_csv(path, sep=r"\s+", header=None,
                     names=["time_ms_file", "au"], engine="python")

    # Drop any rows with missing values, then renumber the rows 0,1,2,...
    # (reset_index). drop=True throws away the old row numbers.
    df = df.dropna().reset_index(drop=True)

    # Compute the true time step between frames, in milliseconds.
    # 8.5 frames/sec -> 1000/8.5 = ~117.6 ms per frame.
    dt = 1000.0 / FPS

    # Build a brand-new, correct time column:
    #   np.arange(len(df)) makes [0, 1, 2, ..., N-1] (one number per row)
    #   multiplying by dt turns frame indices into milliseconds: 0, 117.6, 235.3...
    df["time_ms"] = np.arange(len(df)) * dt

    # Sanity check: compare the file's own time step to ours, and print a note if
    # they differ by more than 0.5%. This confirms the rounding happened and that
    # rebuilding was the right call. len(df) > 1 guards against a 1-row file.
    if len(df) > 1:
        # np.diff gives the differences between consecutive time values;
        # np.median takes the typical (middle) one, robust to odd values.
        file_dt = np.median(np.diff(df["time_ms_file"].values))
        if abs(file_dt - dt) / dt > 0.005:   # abs() = absolute value
            # An f-string (f"...") lets us drop variable values into text with
            # {curly braces}. The :.2f formats a number to 2 decimal places.
            print(f"    [note] file interval {file_dt:.2f} ms vs true {dt:.2f} ms "
                  f"-> rebuilt time axis from FPS={FPS}")

    # Return TWO things: the table and the time step. In Python you can return
    # several values at once; the caller writes  df, dt = load_trace(path).
    return df, dt


# ==========================================================================
# STAGE 3 -- clean the trace: trim the start, fix zero-dropouts
# ==========================================================================
def clean_trace(df):
    """Trim the first TRIM_START_S seconds (if any) and repair 0 values, which
    MUSCLEMOTION sometimes emits as brief dropouts rather than real data."""

    # Keep only rows whose time is at or beyond the trim point. The bit inside
    # the brackets, df["time_ms"] >= TRIM_START_S*1000, produces a column of
    # True/False; df[...] keeps the True rows. .copy() makes an independent copy
    # so later edits don't trigger pandas warnings about editing a "view".
    df = df[df["time_ms"] >= TRIM_START_S * 1000.0].copy()

    # Replace any exact 0 in the 'au' column with NaN ("Not a Number" = missing).
    #   df.loc[ROW_CONDITION, "column"] = value   assigns to just those cells.
    df.loc[df["au"] == 0, "au"] = np.nan

    # Fill those NaNs by linear interpolation (draw a straight line between the
    # good neighbors on each side). limit_direction="both" also handles a gap at
    # the very start or end.
    df["au"] = df["au"].interpolate(limit_direction="both")

    # Renumber rows again after trimming, and return the cleaned table.
    return df.reset_index(drop=True)


# ==========================================================================
# STAGE 4 -- DETREND: remove the slow baseline drift
# --------------------------------------------------------------------------
# This is the single most important step. The raw trace rises slowly over the
# whole recording because the tissue drifts away from the reference frame, not
# because contractions grow. We estimate that slow rise and subtract it so the
# real contractions sit on a flat baseline around zero.
# ==========================================================================
def detrend(df, dt_ms):
    """Estimate the slow drift with a rolling MEDIAN and subtract it, leaving
    the contraction oscillations centered near zero. dt_ms is the ms per frame."""

    # Convert the detrend window from seconds to a number of samples (frames).
    #   seconds * 1000 = ms;  ms / (ms-per-frame) = number of frames.
    #   round(...) then int(...) makes it a whole number.
    n_win = int(round(DETREND_WINDOW_S * 1000.0 / dt_ms))

    # A centered rolling window works best with an ODD length (equal samples on
    # each side of the middle). The "%" is the remainder operator: n_win % 2 == 0
    # is True when n_win is even; if so, add 1 to make it odd.
    if n_win % 2 == 0:
        n_win += 1                # "+= 1" means "add 1 to n_win"

    # Pull the contraction column out as a plain NumPy array of decimals (float).
    au = df["au"].values.astype(float)

    # Compute the rolling MEDIAN: slide a window of length n_win along the trace
    # and take the middle value at each position. The median follows the slow
    # drift while ignoring the brief contraction humps (a median is barely moved
    # by short spikes). This line is written across multiple lines inside
    # parentheses for readability; Python treats it as one statement.
    baseline = (
        pd.Series(au)                                  # wrap array as a pandas Series
        .rolling(n_win, center=True, min_periods=1)    # rolling window, centered
        .median()                                      # take the median in each window
        .values                                        # back to a plain array
    )
    #   center=True     -> the window is centered on each point (not trailing)
    #   min_periods=1   -> still compute a value near the edges with fewer points

    # Store the baseline (for the QC plot) and the detrended signal (raw minus
    # baseline). After this, df["detr"] oscillates around 0.
    df["baseline"] = baseline
    df["detr"] = au - baseline
    return df


# ==========================================================================
# STAGE 5 -- smooth lightly, then detect peaks (= contractions)
# ==========================================================================
def detect_peaks(df, dt_ms):
    """Smooth the detrended trace a little, then find peaks using two rules:
    a minimum spacing between peaks, and a minimum prominence (how far a peak
    stands out from its surroundings)."""

    x = df["detr"].values          # the detrended signal as a plain array

    # --- gentle smoothing ---
    # Convert the smoothing window from seconds to an odd number of samples,
    # with a floor of 5 (the filter needs a window larger than its polynomial).
    w = int(round(SMOOTH_WINDOW_S * 1000.0 / dt_ms))
    if w % 2 == 0:
        w += 1
    w = max(w, 5)                  # max(w, 5): use w, but never smaller than 5

    # Savitzky-Golay smoothing fits a small parabola (polyorder=2) inside each
    # sliding window. It removes jitter while keeping peak shapes truer than a
    # plain moving average would.
    xs = savgol_filter(x, window_length=w, polyorder=2)
    df["smooth"] = xs             # save the smoothed trace for later + the plot

    # --- minimum spacing between peaks, in samples ---
    min_dist = int(round(MIN_PEAK_DISTANCE_S * 1000.0 / dt_ms))
    min_dist = max(min_dist, 1)   # find_peaks requires at least 1

    # --- estimate the noise level, robustly ---
    # We want "typical wiggle size" WITHOUT the big peaks inflating it. The
    # Median Absolute Deviation (MAD) does that: take the median, measure how far
    # each point sits from it, take the median of those distances. Dividing by
    # 0.6745 rescales MAD so it's comparable to a normal standard deviation.
    noise = np.median(np.abs(xs - np.median(xs))) / 0.6745

    # A peak must stand out by at least this much to count.
    prominence = PROMINENCE_SD_MULT * noise

    # find_peaks returns (1) the positions of the peaks and (2) a dictionary of
    # extra info about them ("props"). We keep both.
    peaks, props = find_peaks(xs, distance=min_dist, prominence=prominence)

    # Return the peak positions, the properties, and the noise level (handy to
    # record per file).
    return peaks, props, noise


# ==========================================================================
# STAGE 6 -- turn one recording into summary NUMBERS
# ==========================================================================
def per_minute_counts(df, peaks, n_minutes=10):
    """Count how many detected contractions fall in each 1-minute bin.

    Bins use absolute recording time on the rebuilt axis: minute 1 = 0-60 s,
    minute 2 = 60-120 s, and so on, so the columns line up across files.
    n_minutes=10 is a 'default argument': if the caller doesn't specify, it's 10."""

    # Times (in seconds) of the detected peaks. df["time_ms"].values[peaks] uses
    # the peak positions to pick out just those rows' times.
    pk_times_s = df["time_ms"].values[peaks] / 1000.0

    # Time of the last sample -> where the recording ends (seconds).
    rec_end_s = df["time_ms"].iloc[-1] / 1000.0   # .iloc[-1] = last row's value

    out = {}   # an empty DICTIONARY: stores named values as key -> value pairs

    # Loop over minutes 1..n_minutes. range(1, n+1) yields 1,2,...,n.
    for minute in range(1, n_minutes + 1):
        lo = (minute - 1) * 60.0     # start of this minute in seconds
        hi = minute * 60.0           # end of this minute in seconds
        col = f"contractions_min_{minute}"   # the column name, e.g. "..._min_3"

        if lo >= rec_end_s:
            # This whole minute is past the end of the recording -> we never
            # observed it. Record NaN (missing), NOT 0, so a short recording
            # doesn't look like it had zero contractions in minutes it lacked.
            out[col] = np.nan
        else:
            # Count peaks whose time is in [lo, hi). The parentheses build two
            # True/False arrays; "&" combines them element-by-element (AND);
            # np.sum counts the Trues; int(...) makes it a plain integer.
            out[col] = int(np.sum((pk_times_s >= lo) & (pk_times_s < hi)))

    return out   # e.g. {"contractions_min_1": 5, "contractions_min_2": 4, ...}


def summarize(df, peaks, dt_ms):
    """Collapse a whole recording into a small set of summary numbers, returned
    as a dictionary. Ordered by how trustworthy each is at this low frame rate
    (most trustworthy first)."""

    # Total recording duration in seconds (last time minus first time).
    dur_s = (df["time_ms"].iloc[-1] - df["time_ms"].iloc[0]) / 1000.0
    n = len(peaks)               # number of detected contractions

    # --- FREQUENCY: the main endpoint. Contractions per minute. ---
    # Robust to drift and to the arbitrary units, which is why it's trustworthy.
    freq_per_min = n / (dur_s / 60.0)

    # --- INTERVAL between contractions, and how regular they are ---
    if n >= 2:                    # need at least 2 peaks to have an interval
        pk_times_s = df["time_ms"].values[peaks] / 1000.0
        intervals = np.diff(pk_times_s)      # gaps between consecutive peaks (s)
        mean_ici = np.mean(intervals)        # average gap
        # Coefficient of variation = SD / mean. Lower = more regular rhythm.
        # ddof=1 makes it the sample standard deviation (the usual choice).
        cv_ici = np.std(intervals, ddof=1) / mean_ici
    else:
        # np.nan = "missing/undefined". Can't compute an interval from <2 peaks.
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
        peak_heights = df["smooth"].values[peaks]   # heights on the detrended trace
        med_h = np.median(peak_heights)             # this recording's typical peak
        # Each peak divided by this recording's own median peak, then averaged.
        # The "... if med_h > 0 else np.nan" guards against dividing by zero.
        peak_mean_to_median = np.mean(peak_heights / med_h) if med_h > 0 else np.nan
        raw_amp_median = med_h
    else:
        peak_mean_to_median, raw_amp_median = np.nan, np.nan

    # --- DUTY CYCLE: what fraction of the time is the tissue "active"? ---
    med = np.median(df["smooth"])                       # middle of the signal
    mad = np.median(np.abs(df["smooth"] - med)) / 0.6745  # robust noise level
    thr = med + mad                                     # "active" threshold
    # Fraction of samples above the threshold. np.mean of a True/False array =
    # the proportion of Trues. float(...) makes it a plain decimal.
    duty_cycle = float(np.mean(df["smooth"].values > thr))

    # Build the result DICTIONARY. round(x, k) trims to k decimals for tidiness.
    # The "... if not np.isnan(x) else np.nan" keeps missing values as missing
    # instead of trying to round them.
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

    # Add the 10 per-minute columns to the same dictionary. dict.update(other)
    # merges another dictionary's entries in.
    result.update(per_minute_counts(df, peaks, n_minutes=10))
    return result


# ==========================================================================
# STAGE 7 -- the QC PLOT (make one per file and actually look at them)
# ==========================================================================
def qc_plot(df, peaks, animal, genotype, out_dir):
    """Save a two-panel picture per recording:
       top    = raw trace with the rolling baseline we subtracted
       bottom = detrended+smoothed trace with detected peaks marked
    Reading these tells you if the settings behaved:
       - baseline wiggling with each contraction -> DETREND_WINDOW_S too short
       - peaks missed or noise marked as peaks   -> tune the peak settings"""

    # Make a figure with 2 stacked plots (2 rows, 1 column) that share the x-axis.
    # fig is the whole figure; ax is a list of the two sub-plots, ax[0] and ax[1].
    fig, ax = plt.subplots(2, 1, figsize=(12, 6), sharex=True)

    t = df["time_ms"] / 1000.0    # x-axis in seconds

    # --- top panel: raw data + baseline ---
    ax[0].plot(t, df["au"], lw=0.6, label="raw")            # lw = line width
    ax[0].plot(t, df["baseline"], lw=1.2, label="rolling baseline")
    ax[0].set_ylabel("Contraction (a.u.)")
    ax[0].legend(loc="upper left", fontsize=8)              # the little key
    ax[0].set_title(f"{animal}  ({genotype})  -- raw + baseline")

    # --- bottom panel: detrended trace + detected peaks ---
    ax[1].plot(t, df["smooth"], lw=0.7, label="detrended + smoothed")
    # Plot a red down-triangle ("rv") at each detected peak. We index the time
    # and smoothed arrays with 'peaks' to get just the peak locations.
    ax[1].plot(t.values[peaks], df["smooth"].values[peaks], "rv",
               ms=6, label=f"peaks (n={len(peaks)})")       # ms = marker size
    ax[1].axhline(0, color="k", lw=0.4)                     # horizontal line at 0
    ax[1].set_ylabel("Detrended (a.u.)")
    ax[1].set_xlabel("Time (s)")
    ax[1].legend(loc="upper left", fontsize=8)

    fig.tight_layout()            # tidy spacing so labels don't overlap

    # Build a safe output filename from the animal ID. This keeps letters,
    # numbers, dash, underscore, and dot as-is, and replaces anything else with
    # "_", so the filename can't contain characters that break file systems.
    safe = "".join(c if c.isalnum() or c in "-_." else "_" for c in animal)
    out_path = os.path.join(out_dir, f"QC_{safe}.png")

    fig.savefig(out_path, dpi=110)   # write the image to disk
    plt.close(fig)                   # free the figure's memory (vital in a loop)
    return out_path                  # return where we saved it


# ==========================================================================
# METADATA -- figure out the genotype from the filename
# ==========================================================================
def parse_metadata(fname):
    """Return a dict {'genotype':..., 'animal':...} for one filename.
    Genotype comes from the first matching keyword in GENOTYPE_KEYWORDS. The
    filename (without its extension) is used as the animal identifier."""

    low = fname.lower()          # lowercase so matching ignores capitalization

    genotype = "UNKNOWN"         # default if nothing matches
    # Loop over each (keyword, label) pair. Tuple-unpacking assigns both at once.
    for keyword, label in GENOTYPE_KEYWORDS:
        if keyword in low:       # 'in' tests substring membership
            genotype = label
            break                # stop at the FIRST match (order matters!)

    if genotype == "UNKNOWN":
        # Warn loudly so a mislabeled/odd filename doesn't slip by silently.
        print(f"    [warn] no genotype keyword in {fname} -> UNKNOWN")

    # os.path.splitext splits "name.txt" into ("name", ".txt"); [0] takes "name".
    animal = os.path.splitext(fname)[0]
    return {"genotype": genotype, "animal": animal}


# ==========================================================================
# STAGE 8 -- process ONE file all the way through
# ==========================================================================
def process_file(path, out_dir):
    """Run every stage on a single file and return its summary dictionary."""

    fname = os.path.basename(path)      # just the filename (drop the folders)
    meta = parse_metadata(fname)        # {'genotype':..., 'animal':...}

    # The pipeline, stage by stage. Each line transforms and passes along 'df'.
    df, dt = load_trace(path)           # read + rebuild time axis
    df = clean_trace(df)                # trim + fix dropouts
    df = detrend(df, dt)                # remove slow drift
    peaks, props, noise = detect_peaks(df, dt)   # find contractions
    metrics = summarize(df, peaks, dt)  # compute the summary numbers

    # Add identifying info + a couple of diagnostics to the metrics dictionary.
    metrics.update({
        "animal": meta["animal"],
        "genotype": meta["genotype"],
        "file": fname,
        "dt_ms": round(dt, 2),
        "noise_au": round(noise, 2),
    })

    # Make and save the QC plot for this file.
    qc_path = qc_plot(df, peaks, meta["animal"], meta["genotype"], out_dir)

    # Print a one-line progress summary to the screen.
    print(f"    {meta['animal']} ({meta['genotype']}): "
          f"n={metrics['n_contractions']}, freq={metrics['freq_per_min']}/min, "
          f"CV={metrics['cv_interval']}  ->  {os.path.basename(qc_path)}")

    return metrics


# ==========================================================================
# STAGE 9 -- the DRIVER: find all files, process each, write the summary CSV
# ==========================================================================
def run(input_dir):
    """Top-level routine: discover files, process them all, save a summary table
    and (below) two Prism-ready per-minute tables."""

    files = discover_files(input_dir)

    # If nothing was found, say so clearly and stop. "return None" exits the
    # function early, handing back the value None ("nothing").
    if not files:                        # 'not [empty list]' is True
        print(f"No files matching '{FILE_GLOB}' in {input_dir} "
              f"(after excluding {EXCLUDE_SUBSTRINGS}).")
        return None

    print(f"Found {len(files)} file(s) in {input_dir}:")
    for f in files:                      # print each file we'll process
        print(f"    {os.path.basename(f)}")

    out_dir = input_dir                  # save outputs alongside the inputs
    print("\nProcessing:")

    # Process every file. This list comprehension calls process_file on each path
    # and collects the returned dictionaries into a list called 'rows'.
    rows = [process_file(f, out_dir) for f in files]

    # Turn the list of per-file dictionaries into a table (one row per file).
    summary = pd.DataFrame(rows)

    # Decide the column order: identity, then metrics, then the 10 per-minute
    # columns, then diagnostics. First build the per-minute column names.
    per_min_cols = [f"contractions_min_{m}" for m in range(1, 11)]
    col_order = (
        ["file", "animal", "genotype", "duration_s", "n_contractions",
         "freq_per_min", "mean_interval_s", "cv_interval",
         "peak_mean_to_median", "duty_cycle", "raw_amp_median_au"]
        + per_min_cols                    # "+" joins two lists end to end
        + ["dt_ms", "noise_au"]
    )
    # Reorder columns, but only keep names that actually exist (a safety guard).
    summary = summary[[c for c in col_order if c in summary.columns]]

    # Write the main summary table to a CSV file next to the inputs.
    out_csv = os.path.join(out_dir, "per_recording_summary.csv")
    summary.to_csv(out_csv, index=False)     # index=False: don't write row numbers
    print(f"\nWrote {out_csv}")
    print(summary.to_string(index=False))    # also print it to the screen

    # ---- PRISM-READY PER-MINUTE TABLES -----------------------------------
    # For the per-minute line graph in Prism you want rows = minute (1-10) and
    # columns = animals, with a separate block per genotype. We make one CSV per
    # genotype in that orientation so you can paste each into its Prism group.
    #
    # summary.groupby("genotype") splits the table into sub-tables, one per
    # genotype. The loop variable 'geno' is the genotype label; 'block' is that
    # genotype's sub-table.
    for geno, block in summary.groupby("genotype"):
        # From this genotype's block, take the per-minute columns, using 'animal'
        # as the row label, then .T = transpose (swap rows and columns) so that
        # minutes become rows and animals become columns.
        t = block.set_index("animal")[per_min_cols].T
        t.index = range(1, 11)           # relabel the rows 1..10 (the minutes)
        t.index.name = "minute"
        prism_path = os.path.join(out_dir, f"prism_per_minute_{geno}.csv")
        t.to_csv(prism_path)
        print(f"Wrote {prism_path}  (rows=minute, cols=animals)")

    return summary   # hand the finished table back (useful when run interactively)


# ==========================================================================
# ENTRY POINT -- what actually runs when you launch the script
# --------------------------------------------------------------------------
# Every Python file has a hidden variable __name__. When you RUN this file
# directly (python3 thisfile.py), __name__ equals "__main__", so the block
# below runs. When another file IMPORTS this one, __name__ is the module name
# instead, so the block does NOT run -- letting you reuse the functions without
# triggering a full analysis. This is a standard Python idiom.
# ==========================================================================
if __name__ == "__main__":
    # Set up reading command-line options.
    parser = argparse.ArgumentParser(
        description="Analyze MUSCLEMOTION contraction traces.")
    # Define one option, --input-dir. If omitted, it falls back to INPUT_DIR.
    parser.add_argument("--input-dir", default=INPUT_DIR,
                        help="Directory containing *contraction*.txt files "
                             "(default: current directory).")
    args = parser.parse_args()     # read what the user actually typed
    run(args.input_dir)            # go: run the whole pipeline on that folder


# ==========================================================================
# APPENDIX -- stepping through this interactively (for learning)
# --------------------------------------------------------------------------
# On the server you can run pieces by hand to watch each stage. In an IPython or
# python3 session, load just the definitions (everything above the entry point):
#
#     exec(open("musclemotion-verbose.py").read().split("if __name__")[0])
#
# Then, for one file:
#
#     files = discover_files("/path/to/folder")   # list the files
#     df, dt = load_trace(files[0])               # read the first one
#     df = clean_trace(df)                         # clean it
#     df = detrend(df, dt)                         # remove drift
#     df.head()          # peek at the first rows: columns au/baseline/detr appear
#     peaks, props, noise = detect_peaks(df, dt)   # find peaks
#     peaks              # the frame indices of detected contractions
#     summarize(df, peaks, dt)   # the summary numbers for this one file
#
# Change a CONFIG value at the top, re-run just that stage, and look again. That
# tight loop -- tweak, re-run one stage, inspect -- is the whole point of working
# interactively, and it's how you tune the detrend/peak settings against the QC
# plots rather than against a number you were hoping to see.
# ==========================================================================
