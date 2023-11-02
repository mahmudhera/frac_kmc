import pandas as pd
import matplotlib.pyplot as plt

# Read the data into a pandas DataFrame
data = pd.read_csv("records")

# Get unique scaled values
scaled_values = data['scaled'].unique()

# Create a new figure
fig, axs = plt.subplots(1, len(scaled_values), figsize=(3*len(scaled_values), 3))

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
    if i == 0:
        ax.set_ylabel('Time(s)')
    ax.set_title(f'Scaled = {scaled_value}')

    ax.set_ylim(0, 150)

    # Add a legend
    ax.legend()

# Adjust spacing between subplots
plt.tight_layout()

# Display the figure
#plt.show()
plt.savefig('res.png', format='png', dpi=1200)
