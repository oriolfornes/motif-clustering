import numpy as np
import re

def get_pwms(meme_file):

    # Initialize
    arr = []
    pwms = {}
    r = re.compile(r"\S+")  # Compile a regexp to capture float values

    with open(meme_file) as handle:
        for line in handle:
            if line.startswith("MOTIF"):
                arr.append([line.split(" ")[1], []])
            if len(arr) == 0:
                continue
            try:
                floats = [float(i) for i in r.findall(line)]
                if len(floats) == 4:
                    arr[-1][-1].append(floats)
            except:
                pass

    for k, v in arr:
        pwms.setdefault(k, np.array(v))

    return(pwms)