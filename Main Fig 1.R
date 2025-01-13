library(forestplot)

# Example data
exposures <- c("Seropositivity", "Titers", "Seropositivity", "Titers", "Seropositivity", "Titers","Seropositivity", "Titers")
outcomes <- c("Dementia", "Dementia", "Cognitive Impairment", "Cognitive Impairment", "Alzheimer´s Disease", "Alzheimer´s Disease","Dementia Free Survival","Dementia Free Survival")
ORs <- c(0.93, 0.96, 1.03, 0.98, 1.11, 0.98,1.07,1.04)
lower_CIs <- c(0.66, 0.82, 0.81, 0.88, 0.74, 0.81,0.76,0.88)
upper_CIs <- c(1.31, 1.13, 1.31, 1.10, 1.66, 1.20,1.51,1.23)

data <- data.frame(
  `Cognitive Status` = outcomes,
  `Helicobacter Pylori` = exposures,
  OR = ORs,
  Lower_CI = lower_CIs,
  Upper_CI = upper_CIs
)

# Prepare the table text
tabletext <- rbind(
  c("Cognitive Status", "Helicobacter Pylori", "OR (95% CI)"),
  c("Dementia", "Seropositivity", paste0(format(ORs[1], digits = 3), " (", format(lower_CIs[1], digits = 3), ", ", format(upper_CIs[1], digits = 3), ")")),
  c("", "Titers", paste0(format(ORs[2], digits = 3), " (", format(lower_CIs[2], digits = 3), ", ", format(upper_CIs[2], digits = 3), ")")),
  c("Cognitive Impairment", "Seropositivity", paste0(format(ORs[3], digits = 3), " (", format(lower_CIs[3], digits = 3), ", ", format(upper_CIs[3], digits = 3), ")")),
  c("", "Titers", paste0(format(ORs[4], digits = 3), " (", format(lower_CIs[4], digits = 3), ", ", format(upper_CIs[4], digits = 3), ")")),
  c("Alzheimer's Disease", "Seropositivity", paste0(format(ORs[5], digits = 3), " (", format(lower_CIs[5], digits = 3), ", ", format(upper_CIs[5], digits = 3), ")")),
  c("", "Titers", paste0(format(ORs[6], digits = 3), " (", format(lower_CIs[6], digits = 3), ", ", format(upper_CIs[6], digits = 3), ")")),
  c("Dementia Free Survival", "Seropositivity", paste0(format(ORs[7], digits = 3), " (", format(lower_CIs[7], digits = 3), ", ", format(upper_CIs[7], digits = 3), ")")),
  c("", "Titers", paste0(format(ORs[8], digits = 3), " (", format(lower_CIs[8], digits = 3), ", ", format(upper_CIs[8], digits = 3), ")"))
)

# Define the data structure for forestplot
plot_data <- rbind(
  c(NA, NA, NA), # NA row for header
  cbind(data$OR, data$Lower_CI, data$Upper_CI)
)

# Set the TIFF device
tiff("forestplot_output.tiff", width = 8, height = 6, units = "in", res = 300)

# Create the forest plot
forestplot(
  labeltext = tabletext,
  mean = as.numeric(plot_data[, 1]),
  lower = as.numeric(plot_data[, 2]),
  upper = as.numeric(plot_data[, 3]),
  xlab = "Odds Ratio",
  new_page = TRUE,
  boxsize = 0.15,
  col = fpColors(box = "black", lines = "black", zero = "black"),
  xticks = c(0.5, 1, 1.5, 2),  # Customize x-axis tick marks
  txt_gp = fpTxtGp(label = list(gpar(fontsize = 12, fontfamily = "times")),
                   ticks = gpar(fontsize = 16)),  # Adjust font size of x-axis tick labels
  linewidth = 7,  # Increase line width of the plot
  ci.vertices = TRUE,  # Display vertices for confidence intervals
  ci.vertices.height = 0.05,  # Adjust the height of the CI vertices
  ci.vertices.width = 0.2,  # Adjust the width of the CI vertices
  lwd.ci = 1,  # Line width of the confidence interval lines
  lwd.zero = 2,  # Increased line width of the zero reference line
  zero = 1,  # Position of the zero reference line
  grid = TRUE, 
  boxfill = c("black", "white"),  # Fill color for boxes
  row_margin = unit(0.05, "npc"),  # Further reduce the space between all rows uniformly
  colgap = unit(0.01, "npc")  # Keep the space between columns the same
)

# Turn off the TIFF device
dev.off()


