import os
import json

fname = os.path.join(__path__[0], "../config.json")
if os.path.exists(fname):
    with open(fname, 'r') as f:
        CONFIGS = json.loads(f.read())
else:
    CONFIGS = {}


def test():
    print(
        " _       __     __                             __           __  __                 __      __  "
    )
    print(
        "| |     / /__  / /________  ____ ___  ___     / /_____     / /_/ /_  ___     _____/ /_  __/ /_ "
    )
    print(
        "| | /| / / _ \/ / ___/ __ \/ __ `__ \/ _ \   / __/ __ \   / __/ __ \/ _ \   / ___/ / / / / __ \\"
    )
    print(
        "| |/ |/ /  __/ / /__/ /_/ / / / / / /  __/  / /_/ /_/ /  / /_/ / / /  __/  / /__/ / /_/ / /_/ /"
    )
    print(
        "|__/|__/\___/_/\___/\____/_/ /_/ /_/\___/   \__/\____/   \__/_/ /_/\___/   \___/_/\__,_/_.___/ "
    )
