# =============================================================================
# HASTS 416 — Stochastic Processes | Tutorial 1
# Fully corrected script for A1, A2 and A3
# Plots render in the normal R / RStudio Plots pane

#Mazvita N Joramu
#R228074B
#HASTS

# =============================================================================

# -----------------------------------------------------------------------------
# 0. Packages
# -----------------------------------------------------------------------------
pkgs <- c("markovchain", "igraph", "ggplot2", "reshape2", "scales", "ggrepel", "grid")
invisible(lapply(pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}))

set.seed(2024)

# -----------------------------------------------------------------------------
# 1. Colours and theme
# -----------------------------------------------------------------------------
BLUE   <- "#2166AC"
ORANGE <- "#E08214"
PURPLE <- "#762A83"
RED    <- "#D6604D"
GREEN  <- "#4DAC26"
TEAL   <- "#35978F"
GREY   <- "#5F6368"

base_theme <- theme_classic(base_size = 13) +
  theme(
    plot.title         = element_text(face = "bold", hjust = 0.5, size = 15,
                                      margin = margin(b = 6)),
    plot.subtitle      = element_text(hjust = 0.5, colour = "grey40",
                                      size = 11.3, margin = margin(b = 8)),
    plot.caption       = element_text(colour = "grey45", size = 9.5,
                                      hjust = 0, margin = margin(t = 8)),
    axis.title         = element_text(face = "bold", size = 12),
    axis.text          = element_text(size = 11),
    legend.position    = "bottom",
    legend.title       = element_text(face = "bold"),
    legend.text        = element_text(size = 10.5),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.5),
    panel.grid.major.x = element_line(colour = "grey95", linewidth = 0.4),
    plot.margin        = margin(20, 25, 15, 25)
  )
theme_set(base_theme)

# -----------------------------------------------------------------------------
# 2. Helper functions
# -----------------------------------------------------------------------------
sim_mc <- function(P, n_steps, start = NULL) {
  n <- nrow(P)
  if (is.null(start)) start <- sample(n, 1)
  path <- integer(n_steps + 1)
  path[1] <- start
  for (t in 2:(n_steps + 1)) {
    path[t] <- sample(n, 1, prob = P[path[t - 1], ])
  }
  path
}

uncond_probs <- function(P, pi0, N) {
  k <- ncol(P)
  out <- matrix(NA_real_, N + 1, k)
  out[1, ] <- pi0
  Pn <- diag(k)
  for (t in seq_len(N)) {
    Pn <- Pn %*% P
    out[t + 1, ] <- pi0 %*% Pn
  }
  out
}

# =============================================================================
# A1
# =============================================================================
cat(strrep("=", 78), "\n")
cat("A1 — 5-State Markov Chain\n")
cat(strrep("=", 78), "\n\n")

P_A1 <- matrix(c(
  1.0, 0.0, 0.0, 0.0, 0.0,
  0.5, 0.0, 0.0, 0.0, 0.5,
  0.2, 0.0, 0.0, 0.0, 0.8,
  0.0, 0.0, 1.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 1.0, 0.0
), nrow = 5, byrow = TRUE)

states_A1 <- paste0("S", 1:5)
rownames(P_A1) <- colnames(P_A1) <- states_A1

mc_A1 <- new("markovchain",
             transitionMatrix = P_A1,
             states = states_A1,
             name = "A1: 5-State Markov Chain")

# ──────────────────────────────────────────────────────────────────────────────
#  A1(a): Diagram & Classification
# ──────────────────────────────────────────────────────────────────────────────
cat("─── A1(a): Markov Chain Diagram & Classification ────────────────────────\n\n")

# State classification (analytical):
# ─ S1: ABSORBING  (P[1,1] = 1.0 → self-loop, never leaves)
# ─ S2: TRANSIENT  (exits to S1 or S5; neither can return to S2  →  f₂₂ = 0)
# ─ S3, S4, S5: TRANSIENT with period 3
#   Cycle: S3 → S5 → S4 → S3  (length 3, prob = 0.8×1×1 = 0.8 per visit)
#   Escape: S3 → S1 with probability 0.2 per visit, so eventual absorption.
#
# Periods:
#   S1          → period = 1  (absorbing, aperiodic)
#   S2          → period undefined  (no return path, f₂₂ = 0)
#   S3, S4, S5  → period = 3  (gcd{n ≥ 1 : Pⁿ(i,i) > 0} = gcd{3,6,9,...} = 3)

# Colours: blue = recurrent/absorbing, orange = transient
vcol_A1 <- c(BLUE, ORANGE, ORANGE, ORANGE, ORANGE)

dev.new(width = 7, height = 6, noRStudioGD = TRUE)
draw_mc(P_A1, s1, "A1(a): 5-State Markov Chain", vcol_A1)
legend("bottomright",
       legend = c("Recurrent (Absorbing)", "Transient"),
       fill   = c(BLUE, ORANGE),
       border = "white", bty = "n", cex = 0.88)

cat("Communicating classes:\n")
cc1 <- communicatingClasses(mc_A1)
for (i in seq_along(cc1)) cat("  Class", i, ":", paste(cc1[[i]], collapse = ", "), "\n")

cat("\nRecurrent classes   :", paste(unlist(recurrentClasses(mc_A1)), collapse = ", "), "\n")
cat("Transient states    :", paste(unlist(transientClasses(mc_A1)), collapse = ", "), "\n")
cat("Absorbing states    :", paste(absorbingStates(mc_A1), collapse = ", "), "\n")
cat("Reflecting states   : none\n")

cat("\nPeriods:\n")
cat("  S1          → period = 1  (absorbing / aperiodic)\n")
cat("  S2          → period = undefined (f₂₂ = 0; no return path)\n")
cat("  S3, S4, S5  → period = 3  (shortest return loop: S3→S5→S4→S3)\n")

# -----------------------------------------------------------------------------
# A1(b)
# -----------------------------------------------------------------------------
cat(strrep("-", 78), "\n")
cat("A1(b) — Three simulated trajectories\n")
cat(strrep("-", 78), "\n\n")

N_steps_A1 <- 35
starts_A1 <- sample(1:5, 3, replace = FALSE)

cat("Random starting states selected: ",
    paste(paste0("S", starts_A1), collapse = ", "),
    "\n\n", sep = "")

traj_df_A1 <- do.call(rbind, lapply(seq_along(starts_A1), function(i) {
  data.frame(
    Step  = 0:N_steps_A1,
    State = sim_mc(P_A1, N_steps_A1, start = starts_A1[i]),
    Label = factor(
      paste0("Trajectory ", i, " (start = S", starts_A1[i], ")"),
      levels = paste0("Trajectory ", seq_along(starts_A1), " (start = S", starts_A1, ")")
    )
  )
}))

traj_cols_A1 <- setNames(c(BLUE, RED, GREEN), levels(traj_df_A1$Label))

p_A1b <- ggplot(traj_df_A1, aes(x = Step, y = State, colour = Label, group = Label)) +
  geom_step(linewidth = 1.15) +
  geom_point(size = 2.0, alpha = 0.88) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = BLUE, alpha = 0.30) +
  scale_colour_manual(values = traj_cols_A1, name = NULL) +
  scale_y_continuous(breaks = 1:5, labels = states_A1, expand = expansion(add = 0.35)) +
  scale_x_continuous(breaks = seq(0, N_steps_A1, 5)) +
  labs(
    title = "A1(b): Three Simulated Trajectories",
    subtitle = "Each trajectory eventually enters S1 and remains there",
    x = "Step n",
    y = "State"
  ) +
  guides(colour = guide_legend(override.aes = list(linewidth = 2.5))) +
  theme(legend.key.width = unit(1.8, "cm"))
print(p_A1b)

cat("Comment:\n")
cat("  All three trajectories eventually reach S1 and remain there forever.\n")
cat("  Before absorption, some paths may move through the transient cycle involving S3, S4 and S5.\n\n")

# -----------------------------------------------------------------------------
# A1(c)
# -----------------------------------------------------------------------------
cat(strrep("-", 78), "\n")
cat("A1(c) — Steady-state probabilities and ergodicity\n")
cat(strrep("-", 78), "\n\n")

ss_A1 <- steadyStates(mc_A1)
cat("Stationary distribution π:\n")
print(round(ss_A1, 6))

cat("\nInterpretation:\n")
cat("  π = (1, 0, 0, 0, 0).\n")
cat("  In the long run, all probability mass is concentrated in S1.\n")
cat("  The chain is NOT ergodic because it is not irreducible and contains transient states.\n\n")

# -----------------------------------------------------------------------------
# A1(d)
# -----------------------------------------------------------------------------
cat(strrep("-", 78), "\n")
cat("A1(d) — Unconditional probabilities over time\n")
cat(strrep("-", 78), "\n\n")

N_d_A1 <- 45
pi0_A1 <- rep(1/5, 5)
up_A1 <- uncond_probs(P_A1, pi0_A1, N_d_A1)

df_A1d <- reshape2::melt(
  data.frame(n = 0:N_d_A1, setNames(as.data.frame(up_A1), states_A1)),
  id.vars = "n",
  variable.name = "State",
  value.name = "Probability"
)

end_df_A1 <- subset(df_A1d, n == N_d_A1)
cols_A1 <- c("S1" = BLUE, "S2" = RED, "S3" = GREEN, "S4" = ORANGE, "S5" = PURPLE)

p_A1d <- ggplot(df_A1d, aes(x = n, y = Probability, colour = State)) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = BLUE, linewidth = 0.7, alpha = 0.6) +
  annotate("text", x = N_d_A1 - 2, y = 1.03, label = "π = 1 for S1",
           colour = BLUE, size = 4, hjust = 1, fontface = "italic") +
  ggrepel::geom_label_repel(
    data = end_df_A1,
    aes(label = State),
    size = 3.5,
    fontface = "bold",
    show.legend = FALSE,
    nudge_x = 2.0,
    direction = "y",
    segment.colour = "grey60",
    box.padding = 0.25,
    point.padding = 0.20
  ) +
  scale_colour_manual(values = cols_A1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(-0.01, 1.08)) +
  scale_x_continuous(breaks = seq(0, N_d_A1, 5), expand = expansion(add = c(0.5, 6))) +
  labs(
    title = "A1(d): Unconditional State Probabilities Over Time",
    subtitle = "Initial distribution π0 = (0.2, 0.2, 0.2, 0.2, 0.2)",
    x = "Step n",
    y = "Probability π_n(.)"
  ) +
  theme(legend.position = "none")
print(p_A1d)

cat("Comment on convergence speed:\n")
cat("  Convergence is fairly fast: by about 20 to 25 steps, the unconditional probabilities are already very close to (1, 0, 0, 0, 0).\n\n")

# =============================================================================
# A2
# =============================================================================
cat(strrep("=", 78), "\n")
cat("A2 — 7-State Markov Chain\n")
cat(strrep("=", 78), "\n\n")

P_A2 <- matrix(c(
  0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.4, 0.2, 0.2, 0.2,
  0.0, 0.0, 0.0, 0.0, 0.2, 0.4, 0.4,
  0.3, 0.0, 0.0, 0.1, 0.3, 0.1, 0.2,
  0.0, 0.0, 0.0, 0.2, 0.2, 0.3, 0.3,
  0.0, 0.0, 0.0, 0.5, 0.2, 0.2, 0.1
), nrow = 7, byrow = TRUE)

states_A2 <- paste0("S", 1:7)
rownames(P_A2) <- colnames(P_A2) <- states_A2

mc_A2 <- new("markovchain",
             transitionMatrix = P_A2,
             states = states_A2,
             name = "A2: 7-State Markov Chain")

# ──────────────────────────────────────────────────────────────────────────────
#  A2(a): Diagram
# ──────────────────────────────────────────────────────────────────────────────
cat("─── A2(a): Markov Chain Diagram ─────────────────────────────────────────\n\n")

# Custom layout: {S1,S2} on the left; {S3} at top-centre; {S4–S7} on the right
lo_A2 <- matrix(c(
  -2.0,  0.4,    # S1
  -2.0, -0.4,    # S2
  -0.3,  1.8,    # S3
  1.6,  1.0,    # S4
  1.6,  0.0,    # S5
  1.6, -1.0,    # S6
  -0.3, -1.8     # S7
), ncol = 2, byrow = TRUE)

# Colours:
#   {S1, S2} → BLUE   (recurrent, period 2)
#   {S3}     → PURPLE (transient singleton)
#   {S4–S7}  → ORANGE (transient communicating class)
vcol_A2 <- c(BLUE, BLUE, PURPLE, ORANGE, ORANGE, ORANGE, ORANGE)

dev.new(width = 8, height = 7, noRStudioGD = TRUE)
draw_mc(P_A2, s2, "A2(a): 7-State Markov Chain", vcol_A2,
        layout_fn = lo_A2, curve = 0.28, vsize = 26)
legend("bottomright",
       legend = c("Recurrent {S1,S2}  —  period 2",
                  "Transient singleton {S3}",
                  "Transient class {S4,S5,S6,S7}"),
       fill   = c(BLUE, PURPLE, ORANGE),
       border = "white", bty = "n", cex = 0.82)

# -----------------------------------------------------------------------------
# A2(b)
# -----------------------------------------------------------------------------
cat(strrep("-", 78), "\n")
cat("A2(b) — Classification and period analysis\n")
cat(strrep("-", 78), "\n\n")

cat("Communicating classes:\n")
cc_A2 <- communicatingClasses(mc_A2)
for (i in seq_along(cc_A2)) {
  cat(sprintf("  Class %d : {%s}\n", i, paste(cc_A2[[i]], collapse = ", ")))
}

cat("\nRecurrent classes:\n")
rc_A2 <- recurrentClasses(mc_A2)
for (i in seq_along(rc_A2)) {
  cat(sprintf("  Recurrent class %d : {%s}\n", i, paste(rc_A2[[i]], collapse = ", ")))
}

cat("\nTransient states:\n")
cat("  ", paste(transientStates(mc_A2), collapse = ", "), "\n", sep = "")
cat("\nAbsorbing states:\n  none\n")
cat("\nReflecting states:\n  S5, S6 and S7 are reflecting states because p55 > 0, p66 > 0 and p77 > 0.\n")

cat("\nPeriods:\n")
cat("  S1 : period = 2\n")
cat("  S2 : period = 2\n")
cat("  S3 : period = undefined (no return path)\n")
cat("  S4 : period = 1\n")
cat("  S5 : period = 1\n")
cat("  S6 : period = 1\n")
cat("  S7 : period = 1\n\n")

cat("Interpretation:\n")
cat("  • {S1, S2} is the only closed recurrent class.\n")
cat("  • S3, S4, S5, S6 and S7 are transient.\n")
cat("  • The chain is not ergodic because it is not irreducible and {S1, S2} has period 2.\n\n")

# ──────────────────────────────────────────────────────────────────────────────
#  A2(c): Simulate 2 trajectories
# ──────────────────────────────────────────────────────────────────────────────
cat("\n─── A2(c): Two Simulated Trajectories ───────────────────────────────────\n\n")

N2c     <- 60
starts2 <- sample(7, 2)
cat("Random starting states:", paste(paste0("S", starts2), collapse = ", "), "\n\n")

traj_list2 <- lapply(seq_along(starts2), function(i)
  data.frame(
    Step  = 0:N2c,
    State = sim_mc(P_A2, N2c, start = starts2[i]),
    Traj  = factor(paste0("Trajectory ", i,
                          "  (start = S", starts2[i], ")"))
  )
)
traj_df2 <- do.call(rbind, traj_list2)

p2c <- ggplot(traj_df2, aes(Step, State, colour = Traj, group = Traj)) +
  geom_step(linewidth = 0.85) +
  geom_point(size = 1.8, alpha = 0.78) +
  scale_colour_manual(values = c(BLUE, RED), name = NULL) +
  scale_y_continuous(breaks = 1:7, labels = s2) +
  labs(
    title    = "A2(c): Two Simulated Trajectories — 7-State Chain",
    subtitle = "Each trajectory begins at a randomly selected state",
    x        = "Step  (n)",
    y        = "State",
    caption  = paste0(
      "Starting states: ", paste(paste0("S", starts2), collapse = ", "), ".\n",
      "After a transient phase, trajectories settle into the period-2 oscillation S1 ↔ S2.")
  ) +
  theme(legend.position = "bottom")

print(p2c)

cat("Comment: After visiting the transient states (S3–S7) a finite number of\n")
cat("times, both trajectories are absorbed into the recurrent class {S1, S2}.\n")
cat("Thereafter the chain alternates deterministically: S1→S2→S1→S2→…\n")
cat("The length of the transient phase depends on the starting state.\n")

# -----------------------------------------------------------------------------
# A2(d)
# -----------------------------------------------------------------------------
cat(strrep("-", 78), "\n")
cat("A2(d) — Limiting / stationary probabilities\n")
cat(strrep("-", 78), "\n\n")

ss_A2 <- steadyStates(mc_A2)
cat("Stationary distribution π:\n")
print(round(ss_A2, 6))

cat("\nInterpretation:\n")
cat("  π = (0.5, 0.5, 0, 0, 0, 0, 0).\n")
cat("  The chain spends equal long-run time in S1 and S2, while transient states have limiting probability 0.\n\n")

N_d_A2 <- 70
pi0_A2 <- rep(1/7, 7)
up_A2 <- uncond_probs(P_A2, pi0_A2, N_d_A2)

df_A2d <- reshape2::melt(
  data.frame(n = 0:N_d_A2, setNames(as.data.frame(up_A2), states_A2)),
  id.vars = "n",
  variable.name = "State",
  value.name = "Probability"
)

end_df_A2 <- subset(df_A2d, n == N_d_A2)
cols_A2 <- c("S1" = BLUE, "S2" = RED, "S3" = PURPLE, "S4" = ORANGE, "S5" = GREEN, "S6" = TEAL, "S7" = GREY)

p_A2d <- ggplot(df_A2d, aes(x = n, y = Probability, colour = State)) +
  annotate("rect", xmin = 0, xmax = 28, ymin = -0.005, ymax = 0.18, fill = ORANGE, alpha = 0.06) +
  annotate("text", x = 19, y = 0.13, label = "Transient decay zone",
           colour = ORANGE, size = 3.6, fontface = "italic") +
  geom_line(linewidth = 1.05) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey40", linewidth = 0.65) +
  annotate("text", x = 4, y = 0.515, label = "0.5 reference",
           colour = "grey35", size = 3.4, fontface = "italic", hjust = 0) +
  ggrepel::geom_label_repel(
    data = end_df_A2,
    aes(label = State),
    size = 3.5,
    fontface = "bold",
    show.legend = FALSE,
    nudge_x = 2.2,
    direction = "y",
    segment.colour = "grey60",
    box.padding = 0.25,
    point.padding = 0.15
  ) +
  scale_colour_manual(values = cols_A2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(-0.01, 0.62)) +
  scale_x_continuous(breaks = seq(0, N_d_A2, 10), expand = expansion(add = c(0.5, 8))) +
  labs(
    title = "A2(d): Unconditional State Probabilities Over Time",
    subtitle = "Initial distribution π0 = (1/7, 1/7, 1/7, 1/7, 1/7, 1/7, 1/7)",
    x = "Step n",
    y = "Probability π_n(.)"
  ) +
  theme(legend.position = "none")
print(p_A2d)

cat("Comment:\n")
cat("  S1 and S2 show the period-2 oscillation of the recurrent class.\n")
cat("  The transient states decay toward zero over time.\n")
cat("  The chain is not ergodic.\n\n")

# =============================================================================
# A3
# =============================================================================
cat(strrep("=", 78), "\n")
cat("A3 — Time-dependent traffic Markov chain\n")
cat(strrep("=", 78), "\n\n")

# Corrected first matrix: second row must sum to 1
P1_A3 <- matrix(c(
  0.4, 0.4, 0.2,
  0.3, 0.4, 0.3,
  0.0, 0.1, 0.9
), nrow = 3, byrow = TRUE)

P2_A3 <- matrix(c(
  0.1, 0.5, 0.4,
  0.1, 0.3, 0.6,
  0.0, 0.1, 0.9
), nrow = 3, byrow = TRUE)

traffic_states <- c("Light", "Heavy", "Jammed")
rownames(P1_A3) <- colnames(P1_A3) <- traffic_states
rownames(P2_A3) <- colnames(P2_A3) <- traffic_states

# -----------------------------------------------------------------------------
# A3(a)
# -----------------------------------------------------------------------------
cat(strrep("-", 78), "\n")
cat("A3(a) — Distribution at 6 PM\n")
cat(strrep("-", 78), "\n\n")

initial_A3 <- c(1, 0, 0)
dist_4pm <- initial_A3
for (i in 1:9) dist_4pm <- dist_4pm %*% P1_A3
dist_6pm <- dist_4pm
for (i in 1:6) dist_6pm <- dist_6pm %*% P2_A3

cat(sprintf("P(Light at 6 PM)  = %.6f\n", dist_6pm[1]))
cat(sprintf("P(Heavy at 6 PM)  = %.6f\n", dist_6pm[2]))
cat(sprintf("P(Jammed at 6 PM) = %.6f\n\n", dist_6pm[3]))

# -----------------------------------------------------------------------------
# A3(b)
# -----------------------------------------------------------------------------
cat(strrep("-", 78), "\n")
cat("A3(b) — Simulation verification with 10,000 trajectories\n")
cat(strrep("-", 78), "\n\n")

set.seed(789)
n_sims_A3 <- 10000
final_state_counts <- integer(3)

for (sim in 1:n_sims_A3) {
  current_state <- 1
  for (t in 1:9) current_state <- sample(1:3, 1, prob = P1_A3[current_state, ])
  for (t in 1:6) current_state <- sample(1:3, 1, prob = P2_A3[current_state, ])
  final_state_counts[current_state] <- final_state_counts[current_state] + 1
}

sim_probs_A3 <- final_state_counts / n_sims_A3

cat("Analytical results:\n")
cat(sprintf("  Light  = %.6f\n", dist_6pm[1]))
cat(sprintf("  Heavy  = %.6f\n", dist_6pm[2]))
cat(sprintf("  Jammed = %.6f\n\n", dist_6pm[3]))

cat("Simulation results:\n")
cat(sprintf("  Light  = %.6f\n", sim_probs_A3[1]))
cat(sprintf("  Heavy  = %.6f\n", sim_probs_A3[2]))
cat(sprintf("  Jammed = %.6f\n\n", sim_probs_A3[3]))

cat("Absolute differences:\n")
cat(sprintf("  Light  = %.6f\n", abs(dist_6pm[1] - sim_probs_A3[1])))
cat(sprintf("  Heavy  = %.6f\n", abs(dist_6pm[2] - sim_probs_A3[2])))
cat(sprintf("  Jammed = %.6f\n\n", abs(dist_6pm[3] - sim_probs_A3[3])))

# Evolution plot
dist_matrix_A3 <- matrix(NA, nrow = 16, ncol = 3)
dist_matrix_A3[1, ] <- initial_A3
current <- initial_A3
for (i in 2:10) {
  current <- current %*% P1_A3
  dist_matrix_A3[i, ] <- current
}
for (i in 11:16) {
  current <- current %*% P2_A3
  dist_matrix_A3[i, ] <- current
}

df_A3 <- data.frame(
  Step = 1:16,
  Time = c("1:00","1:20","1:40","2:00","2:20","2:40","3:00","3:20","3:40","4:00",
           "4:20","4:40","5:00","5:20","5:40","6:00"),
  Light = dist_matrix_A3[, 1],
  Heavy = dist_matrix_A3[, 2],
  Jammed = dist_matrix_A3[, 3]
)

df_A3_long <- reshape2::melt(df_A3, id.vars = c("Step", "Time"),
                             variable.name = "State", value.name = "Probability")

p_A3_1 <- ggplot(df_A3_long, aes(x = Step, y = Probability, colour = State, group = State)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  geom_vline(xintercept = 10, linetype = "dashed", colour = BLUE) +
  annotate("text", x = 10.3, y = 0.92, label = "4 PM:\nP1 → P2", colour = BLUE, hjust = 0) +
  scale_colour_manual(values = c("Light" = GREEN, "Heavy" = ORANGE, "Jammed" = RED)) +
  scale_x_continuous(breaks = 1:16, labels = df_A3$Time) +
  labs(
    title = "A3: Traffic State Distribution Over Time",
    x = "Time",
    y = "Probability"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_A3_1)

bar_df_A3 <- data.frame(
  State = traffic_states,
  Analytical = as.numeric(dist_6pm),
  Simulation = as.numeric(sim_probs_A3)
)

bar_df_A3_long <- reshape2::melt(bar_df_A3, id.vars = "State",
                                 variable.name = "Method", value.name = "Probability")

p_A3_2 <- ggplot(bar_df_A3_long, aes(x = State, y = Probability, fill = State)) +
  geom_col(position = "dodge") +
  facet_wrap(~Method) +
  scale_fill_manual(values = c("Light" = GREEN, "Heavy" = ORANGE, "Jammed" = RED)) +
  labs(
    title = "A3: Analytical vs Simulation Results at 6 PM",
    x = "State",
    y = "Probability"
  ) +
  theme(legend.position = "none")

print(p_A3_2)

cat("Conclusion:\n")
cat("  The simulation closely matches the analytical result.\n")
cat("  By 6 PM, the probability of being Jammed is the highest.\n")
cat("  The switch from P1 to P2 after 4 PM increases the likelihood of congestion.\n\n")

cat(strrep("=", 78), "\n")
cat("Tutorial 1 — Complete\n")
cat(strrep("=", 78), "\n")