"""Settings shared by the enumeration programs.

SUPPLEMENTAL_DATA selects the data folder (default: supplemental_data at
the repository root). It supplies the graphs T1, ..., T74 in O1 and the
published covering when E4 reuses it. Quoted passages are read from
programs/paper_text.json. ENUM_WORKERS sets the worker count.
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
PROGRAMS = os.path.dirname(HERE)
DATA = os.environ.get('SUPPLEMENTAL_DATA', os.path.join(os.path.dirname(PROGRAMS), 'supplemental_data'))


def workers():
    return int(os.environ.get('ENUM_WORKERS', max(1, (os.cpu_count() or 2) - 1)))


def catalogue_records():
    """T1, ..., T74 as recorded in item O1: name, graph6, n, m."""
    with open(os.path.join(DATA, 'no_K5_obstructions', 'O1.json')) as f:
        return json.load(f)['graphs']


def paper_text(label):
    with open(os.path.join(PROGRAMS, 'paper_text.json'), encoding='utf-8') as f:
        return json.load(f)[label]
