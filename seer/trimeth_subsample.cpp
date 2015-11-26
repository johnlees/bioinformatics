// Calculates power of seer
// Assumes ~50% of samples have phenotype

// C/C++/C++11 headers
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <iterator>
#include <vector>
#include <unordered_map>
#include <regex>
#include <thread>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

// Armadillo/boost headers
#include <armadillo>

#include <boost/filesystem.hpp>
using namespace boost::filesystem;

// Test size and range set here
const int start_samples = 50;
const int samples_step = 50;
const int end_samples = 3000;

const double repeats = 100;
const unsigned int rand_seed = time(NULL);

const std::string seer_location = "~/installations/pangwas/src/seer";

struct Sample
{
   std::string sample_name;
   int element_present;
};

// Functions
std::vector<int> reservoir_sample(const size_t size, const size_t max_size); // Indices of samples subsampled
std::string cut_struct_mat(const arma::mat& struct_mat, const std::vector<int>& rows); // Extracts only the used rows of the pop_struct matrix
std::string generate_pheno(const std::vector<Sample>& sample_names, const std::vector<int>& kept_indices);
std::string exec(const char* cmd);
int parse_seer(std::string& seer_output);

// Functors
struct seer_hits
{
   seer_hits(const std::string _kmer_input, const std::vector<Sample> _sample_names, const arma::mat _dsm_mat) : kmer_input(_kmer_input), sample_names(_sample_names), dsm_mat(_dsm_mat)
   {
   }

   int operator()(const size_t num_samples) const
   {
      std::vector<int> samples_kept = reservoir_sample(num_samples, sample_names.size());
      std::string pheno_file = generate_pheno(sample_names, samples_kept);
      std::string struct_mat = cut_struct_mat(dsm_mat, samples_kept);

      std::string seer_cmd = seer_location + " -k " + kmer_input + " -p " + pheno_file + " --struct " + struct_mat;
      std::string seer_return = exec(seer_cmd.c_str());
      int num_hits = parse_seer(seer_return);

      // Delete tmp files
      std::remove(pheno_file.c_str());
      std::remove(struct_mat.c_str());

      return num_hits;
   }

   private:
      std::string kmer_input;
      std::vector<Sample> sample_names;
      arma::mat dsm_mat;
};

std::vector<int> reservoir_sample(const size_t size, const size_t max_size)
{
   std::vector<int> sample_indices;

   for (unsigned int i = 1; i <= size; ++i)
   {
      sample_indices.push_back(i);
   }

   for (unsigned int i = size + 1; i<= max_size; ++i)
   {
      unsigned int j = rand() % i + 1;
      if (j <= size)
      {
         sample_indices[j] = i;
      }
   }

   std::sort(sample_indices.begin(), sample_indices.end());
   return sample_indices;
}

std::string generate_pheno(const std::vector<Sample>& sample_names, const std::vector<int>& kept_indices)
{
   path file_name = unique_path();
   std::ofstream pheno_file(file_name.string().c_str());
   if (!pheno_file)
   {
      throw std::runtime_error("Could not write to tmp pheno file");
   }

   for (auto keep_it = kept_indices.begin(); keep_it != kept_indices.end(); ++keep_it)
   {
      // Generate pheno based on OR here
      int pheno = 0;
      if (sample_names[*keep_it - 1].element_present)
      {
         pheno = 1;
      }

      std::string sample_out = sample_names[*keep_it - 1].sample_name;
      pheno_file << sample_out << "\t" << sample_out << "\t" << pheno << "\n";
   }

   return file_name.string();
}

std::string cut_struct_mat(const arma::mat& struct_mat, const std::vector<int>& rows)
{
   arma::uvec keep_rows(rows.size());
   int i = 0;
   for (auto it = rows.begin(); it != rows.end(); ++it)
   {
      keep_rows[i] = *it - 1;
      ++i;
   }

   arma::mat tmp_struct = struct_mat.rows(keep_rows);

   path file_name = unique_path();
   tmp_struct.save(file_name.string(), arma::hdf5_binary);

   return file_name.string();
}

// From http://stackoverflow.com/questions/478898/how-to-execute-a-command-and-get-output-of-command-within-c
// Captures command output
std::string exec(const char* cmd)
{
    FILE* pipe = popen(cmd, "r");
    if (!pipe) return "ERROR";
    char buffer[128];
    std::string result = "";
    while (!feof(pipe)) {
        if (fgets(buffer, 128, pipe) != NULL)
            result += buffer;
    }
    pclose(pipe);
    return result;
}

int parse_seer(std::string& seer_output)
{
   std::stringstream all_output(seer_output);
   std::string line_buf;

   int hits = 0;
   while (std::getline(all_output, line_buf))
   {
      std::stringstream line_in(line_buf);

      std::string kmer;
      double maf, chisq, adjp, beta;
      if (line_in >> kmer >> maf >> chisq >> adjp >> beta)
      {
         ++hits;
      }
   }

   return hits;
}

int main (int argc, char *argv[])
{
   if (argc != 4)
   {
      throw std::runtime_error("Usage is: ./subsample_seer <sample names.txt> <dsm matrix> <kmer file>");
   }

   // Seed random number generator
   srand(rand_seed);

   // Read in samples and element presence
   std::vector<Sample> all_samples;
   std::ifstream sample_file(argv[1]);
   if (!sample_file)
   {
      throw std::runtime_error("Could not open sample file");
   }

   unsigned int present_total = 0;
   while (sample_file)
   {
      Sample sample_read;
      std::string name_buf, present_buf;
      sample_file >> name_buf >> present_buf;

      if (sample_file)
      {
         sample_read.sample_name = name_buf;
         sample_read.element_present = std::stoi(present_buf);

         all_samples.push_back(sample_read);

         if (sample_read.element_present)
         {
            ++present_total;
         }
      }
   }

   // Read in struct matrix
   arma::mat struct_mat;
   std::string dsm_file_name(argv[2]);
   bool mds_loaded = struct_mat.load(dsm_file_name);
   if (!mds_loaded)
   {
      throw std::runtime_error("Could not load mds matrix " + dsm_file_name);
   }

   // kmer file name
   std::string kmer_file_name(argv[3]);

   // Loop over odds ratios, then sample number
   // (const std::vector<std::string> _sample_names, const arma::mat _dsm_mat, const double _OR, const double _MAF, const double _Sr, const size_t _max_samples
   seer_hits run_seer(kmer_file_name, all_samples, struct_mat);
   for (int num_samples = start_samples; num_samples <= end_samples; num_samples += samples_step)
   {
      for (int repeat = 1; repeat <= repeats; ++repeat)
      {
         std::cout << num_samples << "\t" << repeat << "\t" << run_seer(num_samples) << std::endl;
      }
   }


   return 0;
}

