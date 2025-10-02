SPLITS = ["apple", "cherry", "dragonfruit"]

rule all:
    input:
        expand("split/{part}.txt", part=SPLITS)

rule split_file:
    input:
        "complete/split-data.txt"
    output:
        expand("split/{part}.txt", part=SPLITS)
    run:
        # Read full file contents:
        for line in open(input[0]):
            data = line.strip()
            if data in SPLITS:
                with open(f"split/{data}.txt", "w") as out:
                    out.write(data + "\n")