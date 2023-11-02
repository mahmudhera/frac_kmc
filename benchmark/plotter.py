import pandas as pd
import matplotlib.pyplot as plt

# Read the data into a pandas DataFrame
data = pd.read_csv("records")

# Get unique scaled values
scaled_values = data['scaled'].unique()

# Create a new figure
fig, axs = plt.subplots(len(scaled_values), 1, figsize=(4, 2 * len(scaled_values)))

# Iterate through each scaled value and create a subplot
for i, scaled_value in enumerate(scaled_values):
    subset = data[data['scaled'] == scaled_value]

    ax = axs[i]

    # Plot sourmash_time
    ax.plot(subset['kmer_size'], subset['sourmash_time'], label='sourmash_time', marker='o')
    # Plot frackmc_time
    ax.plot(subset['kmer_size'], subset['frackmc_time'], label='frackmc_time', marker='x')

    # Set labels and title
    ax.set_xlabel('kmer_size')
    ax.set_ylabel('Time')
    ax.set_title(f'Scaled = {scaled_value}')

    # Add a legend
    ax.legend()

# Adjust spacing between subplots
plt.tight_layout()

# Display the figure
plt.show()
