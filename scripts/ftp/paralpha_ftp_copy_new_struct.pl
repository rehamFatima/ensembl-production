#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2025] EMBL-European Bioinformatics Institute
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
# limitations under the License

use warnings;
use Getopt::Long;
use FileHandle;
use Scalar::Util qw(looks_like_number);
use File::Spec::Functions qw/catdir/;
use File::Path  qw(make_path remove_tree);
use File::Copy::Recursive qw(dircopy dirmove);
use File::Copy;
use IO::Compress::Bzip2 qw(bzip2 $Bzip2Error);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use File::Remove 'remove';
use File::Copy qw(move);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Getopt::Long qw(:config no_ignore_case);
use Bio::EnsEMBL::Utils::Exception qw(throw);

my $alpha = 'H';

GetOptions(
  'alpha=s' => \$alpha
);

my $log_file = '>/hps/software/users/ensembl/production/ens2020/modenv/production/ftpdumps_new_structure114/copy_overs/uncharted_'.$alpha.'_03.log';
open LOG, $log_file;

make_dirs_for_asm();

close LOG;
exit;


sub assembly_split {

  my $asm = shift();
  my $subdir_len = 3;
  
  my @substrs = split(/\_/,$asm);  
  my @version = split(/\./,$substrs[1]);
  
  # reduce the length of the string by the version
  my @sub_dirs = ();
  
  for (my $begin = 0; $begin <= length($substrs[1]) - ($subdir_len + length($version[1])); $begin += $subdir_len ) {
    # print substr($substrs[1], $begin, $subdir_len)." beg $begin\n";
    my $sub_dir_str = substr($substrs[1], $begin, $subdir_len);
    push (@sub_dirs, $sub_dir_str);
  }

  return $substrs[0], @sub_dirs, $version[1] ;

}


sub make_dirs_for_asm {
  my $og_ftp_path = "/hps/nobackup/flicek/ensembl/production/ensembl_dumps/ftp_mvp/organisms/";
  opendir my $dir, $og_ftp_path or die "Cannot open directory: $!";
  my @specie_files = readdir $dir;
  closedir $dir;
  my $counter = 0;
  my @annotation_sets = ('ensembl','refseq','flybase','community'); # ,'genbank','wormbase'
  my @community_annotation_sets = ('Salk', 'JGI', 'CGD', 'cnag');

  shift (@specie_files); # .
  shift (@specie_files); # ..
  foreach my $specie (@specie_files){
    next unless $specie =~ /$alpha/;
    opendir my $dir, "$og_ftp_path$specie" or print LOG "$specie\tCannot open directory: $!";

    print "\nspecie : $specie\n";
    my @assemblies = readdir $dir;

    $dir = "";
    shift (@assemblies); # .
    shift (@assemblies); # ..
    foreach my $assembly (@assemblies){
      # next unless $assembly =~ /011100615.1/; 
      my @split_asm = assembly_split($assembly);
      $dir = "";
      foreach my $asd (@split_asm){
        $dir = catdir($dir, $asd);
      }
      
      
      opendir my $dir_asm, "$og_ftp_path$specie/$assembly/" or print LOG "$specie\t$assembly\tCannot open directory: $!";
      my @anno_asms = readdir $dir_asm;
      shift (@anno_asms); # .
      shift (@anno_asms); # ..
      my @dates = ();
      my $asm_src = '';
      my $copy_path = '/hps/nobackup/flicek/ensembl/production/ensembl_dumps/ftp_new_structure'; 
      
      if (-d $copy_path.$dir) {
        print "$copy_path$dir filepath exists, skipping $specie $assembly \n";
        next;
      }
      foreach my $an_asm (@anno_asms){
        
        print "anno asm : $an_asm\n";        
        if ( !((grep $_ eq ($an_asm), @community_annotation_sets) || (grep $_ eq ($an_asm), @annotation_sets) || $an_asm eq 'genome' || $an_asm eq 'vep' )){
          print LOG "Unknown annotation source : $an_asm for $specie assembly $assembly\n"; 
        }
        
        if($an_asm eq 'ensembl' || $an_asm eq 'refseq' || $an_asm eq 'flybase' || $an_asm eq 'community'){ # other annotation sources
        
          $asm_src = $an_asm;
          make_path($copy_path.$dir."/$an_asm");
          print "$og_ftp_path$specie/$assembly/$an_asm/\n";
          opendir my $dir_aa, "$og_ftp_path$specie/$assembly/$an_asm/" or print LOG "$specie\t$assembly\tCannot open directory: $!";
          my @genetic_folders = readdir $dir_aa;
          shift (@genetic_folders); # .
          shift (@genetic_folders); # ..
          my $dir_date;
          
          print "gen folder : $genetic_folders[0]\n";
          if ($genetic_folders[0] eq 'geneset' || $genetic_folders[0] eq 'vep' || $genetic_folders[0] eq 'variation' || $genetic_folders[0] eq 'refseq' || $genetic_folders[0] eq 'homology' ){
            opendir $dir_date, "$og_ftp_path$specie/$assembly/$an_asm/$genetic_folders[0]/" or print LOG "$specie\t$assembly\tCannot open directory: $!";
            @dates = readdir $dir_date;
            shift (@dates); # .
            shift (@dates);
          }
          
          foreach my $g_folder (@genetic_folders){
            opendir my $dir_gen, "$og_ftp_path$specie/$assembly/$an_asm/$g_folder/" or print LOG "$specie\t$assembly\tCannot open directory: $!";
            foreach my $date (@dates){
              make_path($copy_path.$dir."/$an_asm/$date/$g_folder");
              if ($g_folder eq 'geneset'){

                dircopy("$og_ftp_path$specie/$assembly/$an_asm/$g_folder/$date",
                        $copy_path.$dir."/$asm_src/$date/$g_folder/");
                gz_to_bgz($copy_path.$dir."/$asm_src/$date/$g_folder/", $specie, $assembly);
              }
              elsif ($g_folder eq 'homology'){
                
                my $homo = 'homology';
                dircopy("$og_ftp_path$specie/$assembly/$an_asm/$g_folder/$date",
                        $copy_path.$dir."/$asm_src/$date/$g_folder/");
                opendir my $dir_hom, $copy_path.$dir."/$asm_src/$date/$g_folder/";
                my @files = readdir ($dir_hom);
                shift (@files); # .
                shift (@files);
                foreach my $file (@files){
                  if ($file =~ $homo){
                    my $homo_date = date_for_homology($assembly);
                    make_path($copy_path.$dir."/$asm_src/$date/$g_folder/$homo_date");
                    move($copy_path.$dir."/$asm_src/$date/$g_folder/$file",
                    $copy_path.$dir."/$asm_src/$date/$g_folder/$homo_date/homology.tsv.gz");
                  }
                }
              }
              elsif ($g_folder eq 'variation'){
                
                my $version_date = date_for_homology($assembly);
                make_path($copy_path.$dir."/$asm_src/$date/$g_folder/$version_date");
                dircopy("$og_ftp_path$specie/$assembly/$an_asm/$g_folder/$date",
                        $copy_path.$dir."/$asm_src/$date/$g_folder/$version_date");
              }
            }            
            closedir $dir_gen;
          }        
          closedir $dir_aa;
        }
        elsif ( $an_asm eq 'genome') {
          print "genome \n";
          foreach my $date (@dates){
            dircopy("$og_ftp_path$specie/$assembly/$an_asm",
                    $copy_path.$dir."/$asm_src/$date/$an_asm/");
            gz_to_bgz($copy_path.$dir."/$asm_src/$date/$an_asm/", $specie, $assembly);
          }          
        }      

        elsif ( $an_asm eq 'vep') {
          opendir my $dir_vep, "$og_ftp_path$specie/$assembly/$an_asm/";
          my @vep_dirs = readdir $dir_vep;
          shift (@vep_dirs); # .
          shift (@vep_dirs);
          foreach my $vep_dir (@vep_dirs){
            if (grep $_ eq ($vep_dir), @annotation_sets){
              # copy files to equivalent anno path e.g e! or refseq
              # get dir/s (usually geneset)
              # again for "date"
              # pop both the above for . n ..
              # copy files to $v_date/$geneset
              dircopy("$og_ftp_path$specie/$assembly/$an_asm/$asm_src/geneset/$dates[0]",
                      $copy_path.$dir."/$asm_src/$dates[0]/geneset/");
              gz_to_bgz($copy_path.$dir."/$asm_src/$dates[0]/geneset/", $specie, $assembly);
              
            }
            elsif ($vep_dir eq 'genome'){
              #for now copy to ensembl/$date/genome
              dircopy("$og_ftp_path$specie/$assembly/$an_asm/genome/",
                      $copy_path.$dir."/$asm_src/$dates[0]/genome/");
              gz_to_bgz($copy_path.$dir."/$asm_src/$dates[0]/genome/", $specie, $assembly);
            }
          }
          close $dir_vep;
        } 
      }
      closedir $dir_asm;
    }
  }
}

sub gz_to_bgz {
  
  # my $h_path = shift ();
  my $dest_path = shift ();
  my $specie = shift ();
  my $assembly = shift ();
  opendir my $dir, $dest_path or print LOG "$specie\t$assembly\t Cannot open directory: $! $dest_path\n" and return;
  my @files = readdir $dir;
  closedir $dir;
  shift (@files);
  shift (@files);


  foreach my $f (@files){
    if (index ($f, "bgz.gzi") != -1 || index ($f, ".embl") != -1) {
      next;
    }
    if (index ($f, ".fai.gz") != -1 || index ($f, ".gzi.gz") != -1 || index ($f, ".csi.gz") != -1 || index ($f, "bgz.gz") != -1 ) {

      print "nonsense extension file : $f\nSubtr: ".$dest_path.substr ($f, 0, -3 )."\n";
      rename $dest_path.$f, $dest_path.substr ($f, 0, -3 );
      next;
    }
    if (index ($f, "tsv") != -1 || index ($f, "xref") != -1 || index ($f, "chain") != -1 || index ($f, "bgz") != -1 || index ($f, "txt") != -1) { # || index ($f, "csi") != -1 || index ($f, "fai") != -1
      next;
    }

    if(index ($f, ".gz") != -1){
      # extract filename minus .gz
      gunzip ("$dest_path$f", $dest_path.(substr ($f, 0, -3 ))) or print LOG "gunzip failed: $GunzipError\n";

      bzip2 "$dest_path".substr ($f, 0, -3 ), "$dest_path".substr ($f, 0, -3 ).".bgz"
            or print LOG "$specie\t$assembly\tbzip2 failed: $Bzip2Error\n";
      if(index ($f, "gff") == -1 && index ($f, "gtf") == -1){
        remove ($dest_path.$f);        
      }
      remove ($dest_path.substr ($f, 0, -3 ));

    }
  }
} 

sub date_for_homology {  
  
  my $assembly = shift();
  my $dbname = 'ensembl_genome_metadata';
  my $dbuser = 'ensro';
  my $dbpass = '';
  my $dbhost = 'mysql-ens-production-1';
  my $dbport = '4721';

  my $prodb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host => $dbhost,
    -port => $dbport,
    -user => $dbuser,
    -dbname => $dbname,
    -pass => $dbpass,
  );

  my $label_query = $prodb->dbc->prepare("SELECT ensembl_release.label \
      FROM genome  \
      JOIN genome_release ON genome.genome_id = genome_release.genome_id \
      JOIN assembly ON genome.assembly_id = assembly.assembly_id \
      JOIN ensembl_release ON genome_release.release_id = ensembl_release.release_id \
      WHERE assembly.accession='$assembly' \
      AND ensembl_release.status='Released' \
      AND ensembl_release.release_type='Partial' ");

  $label_query->execute();
  while (my $label_row= $label_query->fetchrow_arrayref()){
    my ($label) = @$label_row;
    $label =~ s/\-/_/g;
    print "$label\n";
    return $label;
  }

}
