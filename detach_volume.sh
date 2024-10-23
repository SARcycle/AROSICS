#!/bin/bash

# Check if volume ID is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <volume-id>"
    exit 1
fi

# Variables
VOLUME_ID=$1  # Get the volume ID from the first argument
DEVICE_NAME=/dev/nvme1n1  # Device name for unmounting
MOUNT_POINT=/home/ubuntu/assets

# Unmount the volume
echo "Unmounting the volume from $MOUNT_POINT..."
if mountpoint -q "$MOUNT_POINT"; then
	sudo umount "$MOUNT_POINT"
	echo "Volume unmounted successfully."
else
	echo "No mount point found for $MOUNT_POINT."
fi

# Check the volume state
VOLUME_STATE=$(aws ec2 describe-volumes --volume-ids "$VOLUME_ID" --query "Volumes[0].State" --output text)

if [ "$VOLUME_STATE" == "in-use" ]; then
	echo "Detaching the volume..."
	aws ec2 detach-volume --volume-id "$VOLUME_ID" --force

	# Wait for the volume to be detached
	echo "Waiting for the volume to be detached..."
	TIMEOUT=60  # Set timeout in seconds
	elapsed=0   # Initialize elapsed time

	while true; do
		VOLUME_STATE=$(aws ec2 describe-volumes --volume-ids "$VOLUME_ID" --query "Volumes[0].State" --output text)

		if [ "$VOLUME_STATE" == "available" ]; then
			echo "Volume successfully detached."
			break
		fi

		if [ "$elapsed" -ge "$TIMEOUT" ]; then
			echo "Error: Timeout while waiting for volume to detach."
			exit 1
		fi

		sleep 1
		((elapsed++))
	done
else
	echo "Volume is in state '$VOLUME_STATE'. Detachment not required."
fi

