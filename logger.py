import csv
import datetime
import io
import logging


class CsvFormatter(logging.Formatter):
    """
    A custom logging formatter that outputs log records in CSV format.
    """

    def __init__(self):
        """
        Initializes the formatter with a StringIO buffer and a CSV writer.
        """
        super().__init__()
        self.output = io.StringIO()
        self.writer = csv.writer(self.output, quoting=csv.QUOTE_ALL)
        self.fieldnames = ['Level', 'Message', 'InputFileURL', 'DownloadStartTime', 'DownloadSize', 'DownloadDuration',
                           'WrittenFileName', 'ProcessingStartTime', 'ProcessingDuration']

    def format(self, record):
        """
        Formats a log record into a CSV row.

        :param record: The log record to format.
        :return: The formatted CSV row as a string.
        """
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
    """
    Sets up a logger with a CSV formatter and a file handler.

    :param log_file: The file to write log records to.
    :return: The configured logger.
    """
    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger(__name__)
    logging.root.handlers[0].setFormatter(CsvFormatter())
    file_handler = logging.FileHandler(log_file)
    logger.addHandler(file_handler)
    return logger


def download_file(url):
    """
    Simulates downloading a file from a URL.

    :param url: The URL to download from.
    :return: A tuple containing the download start time, size, and duration.
    """
    start_time = datetime.datetime.now()
    download_size = 123456  # Replace with actual download size
    duration = datetime.timedelta(seconds=5)  # Replace with actual duration
    return start_time, download_size, duration


def process_file(input_file):
    """
    Simulates processing a file.

    :param input_file: The file to process.
    :return: A tuple containing the processing start time, duration, and written file name.
    """
    start_time = datetime.datetime.now()
    duration = datetime.timedelta(seconds=10)  # Replace with actual processing duration
    written_file_name = 'output.txt'  # Replace with actual file name
    return start_time, duration, written_file_name