import csv
import datetime
import io
import logging


class CsvFormatter(logging.Formatter):
    def __init__(self):
        super().__init__()
        self.output = io.StringIO()
        self.writer = csv.writer(self.output, quoting=csv.QUOTE_ALL)
        self.fieldnames = ['Level', 'Message', 'InputFileURL', 'DownloadStartTime', 'DownloadSize', 'DownloadDuration',
                           'WrittenFileName', 'ProcessingStartTime', 'ProcessingDuration']

    def format(self, record):
        extra_data = record.__dict__.get('extra', {})
        row = [record.levelname, record.msg]
        for field in self.fieldnames[2:]:
            row.append(extra_data.get(field, ''))
        self.writer.writerow(row)
        data = self.output.getvalue()
        self.output.truncate(0)
        self.output.seek(0)
        return data.strip()


def setup_logger(log_file):
    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger(__name__)
    logging.root.handlers[0].setFormatter(CsvFormatter())
    file_handler = logging.FileHandler(log_file)
    logger.addHandler(file_handler)
    return logger


def download_file(url):
    # Simulate download logic
    start_time = datetime.datetime.now()
    download_size = 123456  # Replace with actual download size
    duration = datetime.timedelta(seconds=5)  # Replace with actual duration
    return start_time, download_size, duration


def process_file(input_file):
    # Simulate processing logic
    start_time = datetime.datetime.now()
    duration = datetime.timedelta(seconds=10)  # Replace with actual processing duration
    written_file_name = 'output.txt'  # Replace with actual file name
    return start_time, duration, written_file_name
