import datetime

outname = "trialxy"

with open("./out/" + outname + ".txt", "w") as text_file:
    text_file.write(str(datetime.datetime.now()))

print(str(datetime.datetime.now())[-1])
