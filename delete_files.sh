#!/bin/bash

# Loop through all subdirectories in the specified path
for dir in /home/ubuntu/assets/S2/*/*/; do
    # Extract the date from the directory name
    folder_date=$(basename "$dir")

    # Check if the folder name is in the YYYYMMDD format
    if [[ $folder_date =~ ^[0-9]{8}$ ]]; then
        # Convert the folder date to a timestamp
        folder_timestamp=$(date -d "$folder_date" +%s)
        # Get the current date timestamp
        current_timestamp=$(date +%s)
        # Calculate the difference in days
        difference=$(( (current_timestamp - folder_timestamp) / 86400 ))

        # If the difference is greater than 7 days, delete the folder
        if [ "$difference" -gt 7 ]; then
            echo "Deleting: $dir"
            rm -rf "$dir"
            #ls "$dir"
        fi
    fi
done

echo "Löschvorgang abgeschlossen."
