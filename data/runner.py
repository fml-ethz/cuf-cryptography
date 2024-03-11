import pathlib
import subprocess

files = pathlib.Path.cwd().glob('./*.fastq.gz')


def cmd(file):
    file_out = str(file.parent) + "/" + file.name.split(".")[0] + ".filtered." + ".".join(file.name.split(".")[1:])
    cmd = f"bbduk.sh in={file} outm=stdout.fq ref=primer.fasta k=21 hdist=3 | bbduk.sh in=stdin.fq out={file_out} ref=primer.fasta ktrim=r interleaved=f k=21 hdist=3 maxlength=21 mininsert=21 minlength=21"

    output = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    log_out = str(file) + ".log"
    with open(log_out, "w") as text_file:
        text_file.write(output.stderr)
    print(output.stderr)



for file in files:
    cmd(file)