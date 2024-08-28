import zipfile

filenames = ["rot/euler.py", "rot/quat.py", 'rot/rotvec.py']

f = zipfile.ZipFile("answer.zip", "w", zipfile.ZIP_DEFLATED)
for filename in filenames:
    f.write(filename)
f.close()