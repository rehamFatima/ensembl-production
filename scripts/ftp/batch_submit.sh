#!/usr/bin/env bash
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#SBATCH -t 00:05:00

#SBATCH --mem=20

#SBATCH -o /homes/reham/CARSScripts/SLURM_Os/%j.out --error /homes/reham/CARSScripts/SLURM_Os/%j.err

for alpha in {A..Z}
do
  sbatch -t 7-00:00:00 --mem=3G --wrap="perl paralpha_ftp_copy_new_struct.pl -alpha $alpha" -o /homes/reham/CARSScripts/SLURM_Os/%j.out --error /homes/reham/CARSScripts/SLURM_Os/%j.err
  echo "Submitted chromosome $alpha"
done
