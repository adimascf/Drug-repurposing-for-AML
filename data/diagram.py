import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Create the figure and axis
fig, ax = plt.subplots(figsize=(10, 6))

# Define positions for boxes and arrows
boxes = {
    "Input Data": (0.1, 0.8),
    "Variant Identification and Validation": (0.4, 0.8),
    "Risk Scoring": (0.4, 0.6),
    "Integration of Non-Genetic Factors": (0.4, 0.4),
    "Statistical or Machine Learning Models": (0.4, 0.2),
    "Risk Stratification": (0.7, 0.4),
    "Personalized Report": (0.7, 0.2),
}

# Draw boxes
for text, pos in boxes.items():
    ax.add_patch(mpatches.FancyBboxPatch(
        (pos[0], pos[1]), 0.3, 0.1,
        boxstyle="round,pad=0.1", edgecolor="black", facecolor="lightblue"))
    ax.text(pos[0] + 0.15, pos[1] + 0.05, text, ha="center", va="center", fontsize=10)

# Draw arrows between boxes
arrows = [
    ("Input Data", "Variant Identification and Validation"),
    ("Variant Identification and Validation", "Risk Scoring"),
    ("Risk Scoring", "Integration of Non-Genetic Factors"),
    ("Integration of Non-Genetic Factors", "Statistical or Machine Learning Models"),
    ("Statistical or Machine Learning Models", "Risk Stratification"),
    ("Statistical or Machine Learning Models", "Personalized Report"),
]

for start, end in arrows:
    start_pos = boxes[start]
    end_pos = boxes[end]
    ax.annotate("", xy=(end_pos[0], end_pos[1] + 0.05), xytext=(start_pos[0] + 0.3, start_pos[1] + 0.05),
                arrowprops=dict(arrowstyle="->", lw=1.5))

# Add title
ax.set_title("Genomic Data Analysis: Trait and Risk Analysis", fontsize=14, pad=20)

# Remove axes for cleaner look
ax.axis("off")

# Show the diagram
plt.tight_layout()
plt.show()
