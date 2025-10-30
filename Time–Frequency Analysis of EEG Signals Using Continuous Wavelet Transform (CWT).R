# ==========================================================
#  EEG Motor Movement/Imagery Dataset (PhysioNet)
#  Subject: S001 | Runs: R01 and R03
# ==========================================================

# ---- Step 1: Install & load edfReader ----
if(!require(edfReader)) install.packages("edfReader")
library(edfReader)

# ---- Step 2: Define base URL and file names ----
base_url <- "https://physionet.org/files/eegmmidb/1.0.0"
subject <- "S001"
#R01 means the subject 1 has baseline,eyes open function and R03 means Subject 1 is doing a task(imagining opening or closing elft or right fist first)
files_to_get <- c("S001R01.edf", "S001R03.edf")

# ---- Step 3: Create data directory ----
dir.create("eeg_data", showWarnings = FALSE)

# ---- Step 4: Download EDF files "S001R01.edf" and "S001R03.edf" in a folder eeg_data ----
for (fn in files_to_get) {
  url <- file.path(base_url, subject, fn)
  dest <- file.path("eeg_data", fn)
  if (!file.exists(dest)) {
    message("Downloading ", fn, " ...")
    download.file(url, destfile = dest, mode = "wb", quiet = FALSE)
  } else {
    message(fn, " already exists locally.")
  }
}

# ---- Step 5: Define a *robust* EDF reader ----
read_edf_file <- function(path) {
  hdr <- readEdfHeader(path)#reads metadata
  sigs <- readEdfSignals(hdr)#reads actual EEG values of each electrode
  
  # Check all possible structures
  if (!is.null(sigs$signal)) {
    siglist <- sigs$signal
    
  } else if (!is.null(sigs$signalMatrix)) {
    mat <- sigs$signalMatrix
    siglist <- lapply(seq_len(ncol(mat)), function(i) mat[, i])
    names(siglist) <- sapply(hdr$signalHeaders, function(h) h$label)
    
  } else if (inherits(sigs, "ebdfSignals") || is.list(sigs)) {
    # Fallback for newer edfReader versions
    siglist <- sigs
  } else {
    stop("Unrecognized structure in 'sigs' — check with str(sigs).")
  }
  
  samplerates <- sapply(hdr$signalHeaders, function(h) h$sampling)
  
  return(list(hdr = hdr, signals = siglist, samplerates = samplerates))
}

# ---- Step 6: Read both EDFs ----
edf1 <- read_edf_file("eeg_data/S001R01.edf")
edf2 <- read_edf_file("eeg_data/S001R03.edf")

# ---- Step 7: Quick check ----
cat("✅ Both files loaded successfully!\n")
cat("File 1 channels:", length(edf1$signals), "\n")
cat("File 2 channels:", length(edf2$signals), "\n")
names(edf1$hdr)
# Handle both list and data.frame structures
if (is.data.frame(edf1$hdr$sHeaders)) {
  ch_names <- edf1$hdr$sHeaders$label
} else if (is.list(edf1$hdr$sHeaders)) {
  ch_names <- sapply(edf1$hdr$sHeaders, function(h) h$label)
} else {
  stop("Unknown sHeaders structure.")
}

print(ch_names[1:10])
names(edf1$hdr$sHeaders)
# Extract sampling rates
if (is.data.frame(edf1$hdr$sHeaders)) {
  fs_all <- edf1$hdr$sHeaders$sRate
} else if (is.list(edf1$hdr$sHeaders)) {
  fs_all <- sapply(edf1$hdr$sHeaders, function(h) h$sRate)
} else {
  stop("Unknown sHeaders structure — cannot extract sampling rate.")
}

# Get unique sampling rates
unique(fs_all)

# Use the first (common) one
fs <- unique(fs_all)[1]#fs is 160
message("Sampling rate (assumed common): ", fs, " Hz")
# Find indices for motor channels
CH1_IDX <- which(ch_names %in% c("C3", "C3..", "C3."))
CH2_IDX <- which(ch_names %in% c("C4", "C4..", "C4."))

# Print results
message("CH1_IDX (C3): ", CH1_IDX)
message("CH2_IDX (C4): ", CH2_IDX)
ch1_idx<-9
ch2_idx<-13
names(edf1)
sig1_baseline <- as.numeric(edf1$signals[[CH1_IDX]]$signal)
sig2_baseline <- as.numeric(edf1$signals[[CH2_IDX]]$signal)
sig1_movement <- as.numeric(edf2$signals[[CH1_IDX]]$signal)
sig2_movement <- as.numeric(edf2$signals[[CH2_IDX]]$signal)
# Trim to same length if necessary
n <- min(length(sig1_baseline), length(sig1_movement))
n
sig1_baseline <- sig1_baseline[1:n]
sig1_movement <- sig1_movement[1:n]
sig2_baseline <- sig2_baseline[1:n]
sig2_movement <- sig2_movement[1:n]
time_vec <- seq(0, by = 1/fs, length.out = n)
head(time_vec, 5)

# Simple pre-filter: remove DC (mean) 
sig1_baseline <- sig1_baseline - mean(sig1_baseline)
sig1_movement <- sig1_movement - mean(sig1_movement)
sig2_baseline <- sig2_baseline - mean(sig2_baseline)
sig2_movement <- sig2_movement - mean(sig2_movement)
#Continuous Wavelet Transform (CWT) & scalogram (WaveletComp)
# WaveletComp expects a data.frame with a time series column; dt = 1/fs
library(WaveletComp)
# Build data.frame for WaveletComp
df_sig1_baseline <- data.frame(t = time_vec, x = sig1_baseline)
wt_sig1_baseline <- analyze.wavelet(df_sig1_baseline, "x", loess.span = 0,
                                    dt = 1/fs, dj = 1/20, lowerPeriod = 1/(fs/2),
                                    upperPeriod = 1, make.pval = FALSE)
# Plot scalogram
wt.image(wt_sig1_baseline, main = paste("Scalogram - Subject","ch", ch1_idx, "baseline"))
df_sig1_movement <- data.frame(t = time_vec, x = sig1_movement)
wt_sig1_movement <- analyze.wavelet(df_sig1_movement, "x", loess.span = 0,
                                    dt = 1/fs, dj = 1/20, lowerPeriod = 1/(fs/2),
                                    upperPeriod = 1, make.pval = FALSE)
wt.image(wt_sig1_movement, main = paste("Scalogram - Subject","ch", ch1_idx, "movement"))
#Band power (delta/theta/alpha/beta) using bandpass + windowed RMS
library(signal)
bandpass_rms <- function(signal, fs, low, high, win_size_sec = 1, step_sec = 0.5){
  # design Butterworth bandpass (4th order)
  W <- butter(4, c(low, high)/(fs/2), type = "pass")
  filt <- filtfilt(W, signal)
  # window params
  win_pts <- round(win_size_sec * fs)#win_pts means how many data points are in each window
  step_pts <- round(step_sec * fs)#how many points to move the window each time
  starts <- seq(1, length(filt) - win_pts + 1, by = step_pts)# a seq of start indices for each window 
  pow <- sapply(starts, function(s){
    seg <- filt[s:(s+win_pts-1)]
    mean(seg^2)  # power estimate (mean-square)
  })
  return(list(power = pow, starts = starts))
}

# frequency bands (Hz)
bands <- list(delta = c(1,4), theta = c(4,8), alpha = c(8,13), beta = c(13,30))
win_size <- 2    # 2-second windows
step_size <- 1   # 1-second hop

# compute band powers for ch1 baseline and movement
bp_baseline <- lapply(bands, function(b) bandpass_rms(sig1_baseline, fs, b[1], b[2], win_size, step_size)$power)
bp_movement <- lapply(bands, function(b) bandpass_rms(sig1_movement, fs, b[1], b[2], win_size, step_size)$power)
# convert to data.frame for plotting/comparison
nwin <- length(bp_baseline[[1]])
df_bp <- data.frame(window = seq_len(nwin),
                    baseline_alpha = bp_baseline$alpha,
                    movement_alpha = bp_movement$alpha)
# quick plot
library(ggplot2)
ggplot(df_bp, aes(x = window)) +
  geom_line(aes(y = baseline_alpha, color = "baseline")) +
  geom_line(aes(y = movement_alpha, color = "movement")) +
  labs(title = "Alpha band power (example) — channel 1", y = "Mean-square power")
#Wavelet coherence (time-frequency connectivity) between two channels
# Build a data frame with two series
# WaveletComp analyze.coherency
df_two_baseline<-data.frame(t=time_vec,ch1=sig1_baseline,ch2=sig2_baseline)
wc_baseline <- analyze.coherency(df_two_baseline, my.pair = c("ch1","ch2"),
                                 dt = 1/fs, make.pval = FALSE, loess.span = 0)
# Plot coherence image
wc.image(wc_baseline, main = "Wavelet coherence - baseline ch1 vs ch2")

# Movement
df_two_movement <- data.frame(t = time_vec, ch1 = sig1_movement, ch2 = sig2_movement)
wc_movement <- analyze.coherency(df_two_movement, my.pair = c("ch1","ch2"),
                                 dt = 1/fs, make.pval = FALSE, loess.span = 0)
wc.image(wc_movement, main = "Wavelet coherence - movement ch1 vs ch2")
#Phase Locking Value (PLV) across windows (time-resolved)
inst_phase <- function(x) {
  N <- length(x)
  Xf <- fft(x)
  h <- rep(0, N)
  if (N %% 2 == 0) {
    h[1] <- 1
    h[N/2 + 1] <- 1
    h[2:(N/2)] <- 2
  } else {
    h[1] <- 1
    h[2:((N + 1)/2)] <- 2
  }
  x_hilbert <- fft(Xf * h, inverse = TRUE) / N
  Arg(x_hilbert)
}
phase1_baseline <- inst_phase(sig1_baseline)
phase2_baseline <- inst_phase(sig2_baseline)
phase1_movement <- inst_phase(sig1_movement)
phase2_movement <- inst_phase(sig2_movement)
install.packages("ggplot2")
library(gglot2)
# compute PLV per sliding window
compute_plv_windows <- function(phase1, phase2, fs, win_sec=2, step_sec=1){
  win_pts <- round(win_sec * fs)
  step_pts <- round(step_sec * fs)
  starts <- seq(1, length(phase1)-win_pts+1, by = step_pts)
  plv <- sapply(starts, function(s){
    p1 <- phase1[s:(s+win_pts-1)]
    p2 <- phase2[s:(s+win_pts-1)]
    vals <- exp(1i * (p1 - p2))
    abs(mean(vals))
  })
  return(plv)
}
# ---------------------------------------------------------------
# Compute PLV for baseline and movement phases
# ---------------------------------------------------------------

# Baseline phase PLV
plv_baseline <- compute_plv_windows(
  phase1 = phase1_baseline,
  phase2 = phase2_baseline,
  fs = fs,
  win_sec = win_size,
  step_sec = step_size
)

# Movement phase PLV
plv_movement <- compute_plv_windows(
  phase1 = phase1_movement,
  phase2 = phase2_movement,
  fs = fs,
  win_sec = win_size,
  step_sec = step_size
)


# quick plot
df_plv <- data.frame(window = seq_along(plv_baseline),
                     baseline = plv_baseline,
                     movement = plv_movement)
library(tidyr)
df_plv_l <- pivot_longer(df_plv, cols = c("baseline","movement"), names_to = "cond", values_to = "plv")
ggplot(df_plv_l, aes(x = window, y = plv, color = cond)) + geom_line() + labs(title = "PLV across windows")
#Compare connectivity across conditions (paired t-test + simple permutation test)
# We'll compare mean PLV across windows (paired by window index)
# trim to equal length
L <- min(length(plv_baseline), length(plv_movement))
plv_b <- plv_baseline[1:L]; plv_m <- plv_movement[1:L]

# paired t-test
tt <- t.test(plv_b, plv_m, paired = TRUE)
print(tt)

# permutation test (paired): shuffle condition labels across windows
perm_test_paired <- function(x, y, nperm = 5000){
  obs_diff <- mean(y - x)
  n <- length(x)
  combined <- cbind(x,y)
  diffs <- numeric(nperm)
  for(i in seq_len(nperm)){
    # flip sign per pair with 50% prob => equivalent to permuting within pairs
    flips <- sample(c(1,-1), n, replace = TRUE)
    diffs[i] <- mean(flips * (y - x))
  }
  p_val <- mean(abs(diffs) >= abs(obs_diff))
  return(list(obs = obs_diff, p = p_val, null = diffs))
}
perm_res <- perm_test_paired(plv_b, plv_m, nperm = 2000)
print(perm_res)
#Summary outputs to save
saveRDS(list(wt_sig1_baseline = wt_sig1_baseline,
             wt_sig1_movement = wt_sig1_movement,
             wc_baseline = wc_baseline,
             plv_baseline = plv_baseline,
             plv_movement = plv_movement,
             bp_baseline = bp_baseline,
             bp_movement = bp_movement),
        file = "eeg_analysis_results_S001.rds")

message("Done. Results saved to eeg_analysis_results_S001.rds")