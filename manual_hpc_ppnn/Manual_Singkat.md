# Login ke HPC PPNN

Alamat IP publik dari headnode adalah 167.205.6.67.

Anda dapat mengaksesnya melalui SSH jika telah mendapatkan username
dari admin HPC PPNN

```bash
ssh username@167.205.6.67
```

atau

```bash
ssh -l username 167.205.6.67
```

Silakan menghubungi admin jika menemui kesulitan untuk login.

# Submit job

Sistem HPC PPNN menggunakan PBS Pro sebagai job manager.

Contoh script untuk submit job untuk program yang menggunakan MPI.

```bash
#!/bin/bash
#PBS -q normal
#PBS -l nodes=1:ppn=2

cd $PBS_O_WORKDIR

BIN="/app/bin/vasp_5.2_intelmpi_ifort"
mpirun -n 2 $BIN > LOG1
```

Contoh job Gaussian09:

```bash
#!/bin/bash
#PBS -q normal
#PBS -l select=1:ncpus=2

cd $PBS_O_WORKDIR

export g09root=/app
export GAUSS_EXEDIR=$g09root/g09
export PATH=$GAUSS_EXEDIR:$PATH

# NOTES: Please make sure that ncpus is the same as number of processors
# that you requested in the G09 input files.
/app/g09/g09 test0001.com
```