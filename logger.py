import datetime
import hashlib
import json
import os.path


class log:
    def parse_logfile(file_path):
        if os.path.exists(file_path):
            with open(file_path) as f:
                data = json.load(f)
        else:
            data = []
        return data

    def write_logfile(self, file_path):
        pass

    def get_log_entry(self, scene_date, rel_orbit):
        if self.check_log_entry_present(scene_date, rel_orbit) != -1:
            pass

    def set_log_entry(self, scene_date, rel_orbit, property, value):
        pass

    def check_log_entry_present(self, scene_date, rel_orbit, key='hash', value=None):
        if value is None:
            value = hashlib.md5(f'{scene_date.strftime("%Y%m%d")}_{rel_orbit}'.encode('utf-8')).hexdigest()
        for idx, entry in enumerate(self):
            if key in entry and entry[key] == value:
                return idx
        return -1


data = log.parse_logfile('/mnt/d/SATROMO/AROSICS_Coregistration/AROSICS/assets/base_data/log.json')
print(log.check_log_entry_present(data, datetime.date(2022, 3, 1), rel_orbit=65))
