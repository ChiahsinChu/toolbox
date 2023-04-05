import os
import json
import logging

fname = os.path.join(__path__[0], "../config.json")
if os.path.exists(fname):
    with open(fname, 'r') as f:
        CONFIGS = json.loads(f.read())
else:
    CONFIGS = {}

logging.basicConfig(filename="dp_qeq.log",
                    level=logging.INFO,
                    format='%(asctime)s - %(message)s',
                    datefmt='%d-%m-%y %H:%M:%S')
