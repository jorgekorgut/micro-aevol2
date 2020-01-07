//
// Created by arrouan on 01/10/18.
//

#include "Dna.h"
#include "ExpManager.h"

Dna::Dna(const Dna& clone) : seq_(clone.seq_) {
}

Dna::Dna(int length, Threefry::Gen& rng) : seq_(length) {
  // Generate a random genome
  for (int32_t i = 0; i < length; i++) {
    seq_[i] = '0' + rng.random(NB_BASE);
  }
}

Dna::Dna(char* genome, int length) : seq_(length) {
  strcpy(seq_.data(), genome);
}

Dna::Dna(int length) : seq_(length) {
}

int Dna::length() const {
  return seq_.size();
}

void Dna::save(gzFile backup_file) {
    int dna_length = length();
    gzwrite(backup_file, &dna_length, sizeof(dna_length));
    gzwrite(backup_file, seq_.data(), dna_length * sizeof(seq_[0]));
}

void Dna::load(gzFile backup_file) {
    int dna_length;
    gzread(backup_file, &dna_length, sizeof(dna_length));

    char tmp_seq[dna_length];
    gzread(backup_file, tmp_seq, dna_length * sizeof(tmp_seq[0]));

    seq_ = std::vector<char>(tmp_seq, tmp_seq+dna_length);
}

void Dna::set(int pos, char c) {
  seq_[pos] = c;
}

/**
 * Remove the DNA inbetween pos_1 and pos_2
 *
 * @param pos_1
 * @param pos_2
 */
void Dna::remove(int pos_1, int pos_2) {
  assert(pos_1 >= 0 && pos_2 >= pos_1 && pos_2 <= seq_.size());
  seq_.erase(seq_.begin() + pos_1, seq_.begin() + pos_2);
}

/**
 * Insert a sequence of a given length at a given position into the DNA of the Organism
 *
 * @param pos : where to insert the sequence
 * @param seq : the sequence itself
 * @param seq_length : the size of the sequence
 */
void Dna::insert(int pos, std::vector<char> seq) {
// Insert sequence 'seq' at position 'pos'
  assert(pos >= 0 && pos < seq_.size());

  seq_.insert(seq_.begin() + pos, seq.begin(), seq.end());
}

/**
 * Insert a sequence of a given length at a given position into the DNA of the Organism
 *
 * @param pos : where to insert the sequence
 * @param seq : the sequence itself
 * @param seq_length : the size of the sequence
 */
void Dna::insert(int pos, Dna* seq) {
// Insert sequence 'seq' at position 'pos'
  assert(pos >= 0 && pos < seq_.size());

  seq_.insert(seq_.begin() + pos, seq->seq_.begin(), seq->seq_.end());
}

void Dna::do_switch(int pos) {
  if (seq_[pos] == '0') seq_[pos] = '1';
  else seq_[pos] = '0';
}

void Dna::do_duplication(int pos_1, int pos_2, int pos_3) {
  // Duplicate segment [pos_1; pos_2[ and insert the duplicate before pos_3
  char* duplicate_segment = NULL;

  int32_t seg_length;

  if (pos_1 < pos_2) {
    //
    //       pos_1         pos_2                   -> 0-
    //         |             |                   -       -
    // 0--------------------------------->      -         -
    //         ===============                  -         - pos_1
    //           tmp (copy)                      -       -
    //                                             -----      |
    //                                             pos_2    <-'
    //
    std::vector<char> seq_dupl =
        std::vector<char>(seq_.begin() + pos_1, seq_.begin() + pos_2);

    insert(pos_3, seq_dupl);
  } else { // if (pos_1 >= pos_2)
    // The segment to duplicate includes the origin of replication.
    // The copying process will be done in two steps.
    //
    //                                            ,->
    //    pos_2                 pos_1            |      -> 0-
    //      |                     |                   -       - pos_2
    // 0--------------------------------->     pos_1 -         -
    // ======                     =======            -         -
    //  tmp2                        tmp1              -       -
    //                                                  -----
    //
    //
    std::vector<char>
        seq_dupl = std::vector<char>(seq_.begin() + pos_1, seq_.end());
    seq_dupl.insert(seq_dupl.end(), seq_.begin(), seq_.begin() + pos_2);

    insert(pos_3, seq_dupl);
  }
}

int Dna::promoter_at(int pos) {
  int prom_dist[22];

  for (int motif_id = 0; motif_id < 22; motif_id++) {
    // Searching for the promoter
    prom_dist[motif_id] =
        PROM_SEQ[motif_id] ==
        seq_[
            pos + motif_id >= seq_.size() ? pos +
                                            motif_id -
                                            seq_.size()
                                          : pos +
                                            motif_id]
        ? 0
        : 1;

  }


  // Computing if a promoter exists at that position
  int dist_lead = prom_dist[0] +
                  prom_dist[1] +
                  prom_dist[2] +
                  prom_dist[3] +
                  prom_dist[4] +
                  prom_dist[5] +
                  prom_dist[6] +
                  prom_dist[7] +
                  prom_dist[8] +
                  prom_dist[9] +
                  prom_dist[10] +
                  prom_dist[11] +
                  prom_dist[12] +
                  prom_dist[13] +
                  prom_dist[14] +
                  prom_dist[15] +
                  prom_dist[16] +
                  prom_dist[17] +
                  prom_dist[18] +
                  prom_dist[19] +
                  prom_dist[20] +
                  prom_dist[21];

  return dist_lead;
}

int Dna::terminator_at(int pos) {
  int term_dist[4];
  for (int motif_id = 0; motif_id < 4; motif_id++) {

    // Search for the terminators
    term_dist[motif_id] =
        seq_[
            pos + motif_id >= seq_.size() ? pos +
                                            motif_id -
                                            seq_.size() :
            pos + motif_id] !=
        seq_[
            pos - motif_id + 10 >= seq_.size() ?
            pos - motif_id + 10 - seq_.size() :
            pos -
            motif_id +
            10] ? 1
                : 0;
  }
  int dist_term_lead = term_dist[0] +
                       term_dist[1] +
                       term_dist[2] +
                       term_dist[3];

  return dist_term_lead;
}

bool Dna::shine_dal_start(int pos) {
  bool start = false;
  int t_pos, k_t;

  for (int k = 0; k < 9; k++) {
    k_t = k >= 6 ? k + 4 : k;
    t_pos = pos + k_t >= seq_.size() ? pos + k_t -
                                       seq_.size()
                                     : pos + k_t;

    if (seq_[t_pos] ==
        SHINE_DAL_SEQ[k]) {
      start = true;
    } else {
      start = false;
      break;
    }
  }

  return start;
}

bool Dna::protein_stop(int pos) {
  bool is_protein;
  int t_k;

  for (int k = 0; k < 3; k++) {
    t_k = pos + k >= seq_.size() ?
          pos - seq_.size() + k :
          pos + k;

    if (seq_[t_k] ==
        PROTEIN_END[k]) {
      is_protein = true;
    } else {
      is_protein = false;
      break;
    }
  }

  return is_protein;
}

int Dna::codon_at(int pos) {
  int value = 0;

  int t_pos;

  for (int i = 0; i < 3; i++) {
    t_pos =
        pos + i >= seq_.size() ? pos + i -
                                 seq_.size()
                               : pos + i;
    if (seq_[t_pos] ==
        '1')
      value += 1 << (CODON_SIZE - i - 1);
  }

  return value;
}