import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Step 1: Create example data for 3 pairs
data = {
    'Disease': ['Breast Cancer', 'Breast Cancer', 'Inflammatory Bowel Disease', 'Inflammatory Bowel Disease',
                 'Type 2 Diabetes', 'Type 2 Diabetes'],
    'Method': ['GREAT', 'E-MAGO'] * 3,
    'Percentage': [9/150*100, 12/150*100, 12/250*100, 19/110*100, 25/275*100, 25/275*100]  # Example values
}

df = pd.DataFrame(data)

# Step 2: Plotting
plt.figure(figsize=(8, 6))
sns.barplot(data=df, x='Disease', y='Percentage', hue='Method', palette='Set2')

# Step 3: Customize
plt.title('Overlapping genes GREAT vs E-MAGO')
plt.ylabel('Percentage Annotated Genes Overlapping with GC genes (%)')
plt.xlabel('Disease')
plt.ylim(0, 20)
plt.legend(title='Annotation Method')
plt.tight_layout()

# Step 4: Show plot
plt.show()
