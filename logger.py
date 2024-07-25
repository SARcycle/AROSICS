import csv
import hashlib
from datetime import datetime


def write_progress(csv_file, date, orbit_number, column_name, value):
    """
    Write progress to a specific column in a CSV file
    """
    date_str = date.strftime('%Y-%m-%d')
    row_id = hashlib.sha256(f"{date_str}_{orbit_number}".encode()).hexdigest()
    with open(csv_file, 'r+', newline='') as file:
        reader = csv.reader(file)
        writer = csv.writer(file)
        rows = list(reader)
        header_row = rows[0]
        if column_name not in header_row:
            header_row.append(column_name)
            rows[0] = header_row
        for row in rows:
            if row[0] == row_id:
                row[header_row.index(column_name)] = value
                break
        else:
            # If the row doesn't exist, create a new one
            new_row = [row_id] + ["" for _ in range(len(header_row) - 1)]
            new_row[header_row.index(column_name)] = value
            rows.append(new_row)
        writer.writerows(rows)


def read_progress(csv_file, date, orbit_number, column_name):
    """
    Read progress from a specific column in a CSV file
    """
    date_str = date.strftime('%Y-%m-%d')
    row_id = hashlib.sha256(f"{date_str}_{orbit_number}".encode()).hexdigest()
    with open(csv_file, 'r', newline='') as file:
        reader = csv.reader(file)
        header_row = next(reader)
        if column_name not in header_row:
            return None
        for row in reader:
            if row[0] == row_id:
                return row[header_row.index(column_name)]
        return None


# Example usage:
csv_file = '/mnt/d/SATROMO/AROSICS_Coregistration/AROSICS/assets/base_data/log.csv'
date = datetime(2023, 3, 15)
orbit_number = 123
column_name = 'status'

# Write progress
write_progress(csv_file, date, orbit_number, column_name, 'in_progress')

# Read progress
progress = read_progress(csv_file, date, orbit_number, column_name)
print(f"Progress: {progress}")
