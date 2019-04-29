#!/bin/bash

#SBATCH --job-name=smbh		#Nombre del job
#SBATCH -p short			#Cola a usar, Default=short (Ver colas y límites en /hpcfs/shared/README/partitions.txt)
#SBATCH -N 1				#Nodos requeridos, Default=1
#SBATCH -n 1				#Tasks paralelos, recomendado para MPI, Default=1
#SBATCH --cpus-per-task=16		#Cores requeridos por task, recomendado para multi-thread, Default=1
#SBATCH --mem-per-cpu=2000		#Memoria en Mb por CPU, Default=2048
#SBATCH --time=48:00:00			#Tiempo máximo de corrida, Default=2 horas
#SBATCH --mail-user=js.barbosa10@uniandes.edu.co
#SBATCH --mail-type=ALL			
#SBATCH -o smbh.o%j			#Nombre de archivo de salida

module load anaconda/python3
python initial_conditions.py
