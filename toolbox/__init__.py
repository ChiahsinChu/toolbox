import os
import json

fname = os.path.join(__path__[0], "../config.json")
if os.path.exists(fname):
    with open(fname, 'r') as f:
        CONFIGS = json.loads(f.read())
else:
    CONFIGS = {}
