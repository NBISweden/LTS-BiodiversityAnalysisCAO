from datetime import datetime
from math import floor, ceil
from sys import stderr


def title2log(title, logfile, llen=90, also_stderr=True):
    text_insert = "{title} started at : {time}".format(title=title, time=datetime.now())
    prefix = "=" * floor((llen - 2 - len(text_insert)) / 2) + " "
    sufffix = " " + "=" * ceil((llen - 2 - len(text_insert)) / 2)
    text = prefix + text_insert + sufffix
    with open(logfile, "a") as handle:
        handle.writelines("\n\n" + text + "\n\n")
    if also_stderr:
        print(text, file=stderr, flush=True)


def freetxt_line(text, logfile, llen=90, also_stderr=True):
    text_insert = text
    prefix = "=" * floor((llen - 2 - len(text_insert)) / 2) + " "
    sufffix = " " + "=" * ceil((llen - 2 - len(text_insert)) / 2)
    text = prefix + text_insert + sufffix
    with open(logfile, "a") as handle:
        handle.writelines("\n\n" + text + "\n\n")

    if also_stderr:
        print(text, file=stderr, flush=True)
