#pragma once


__device__
void dna_gpu_copy(char** dna, char** dna_next_gen, size_t dna_size, int indiv_id) {
    for (int i = 0; i < dna_size;i++)
        dna_next_gen[indiv_id][i] = dna[indiv_id][i];
}

__device__
void dna_gpu_set(char** dna, int indiv_id, int pos, char c) {
    dna[indiv_id][pos] = c;
}

__device__
void dna_gpu_remove(char** dna, char** dna_next_gen, int indiv_id, size_t dna_size, int pos_1, int pos_2) {
    //assert(pos_1 >= 0 && pos_2 >= pos_1 && pos_2 <= dna_size);

    for(int i = 0; i < pos_1; i++)
        dna_next_gen[indiv_id][i] = dna[indiv_id][i];

    int pos=pos_1;
    for (int i=pos_1; i < dna_size-pos_2;i++,pos++) {
        dna_next_gen[indiv_id][pos] = dna[indiv_id][i];
    }
}


__device__
void dna_gpu_insert(char** dna, char** dna_next_gen, int indiv_id, size_t dna_size, int pos, char* seq, int nb_insert) {
// Insert sequence 'seq' at position 'pos'
    //assert(pos >= 0 && pos < dna_size);

    for (int i = 0; i < pos; i++)
        dna_next_gen[indiv_id][i] = dna[indiv_id][i];

    for (int i = pos, j=0; i < pos+nb_insert; j++,i++)
        dna_next_gen[indiv_id][i] = seq[j];

    for (int i=pos+nb_insert, j=pos; j < dna_size; i++,j++)
        dna_next_gen[indiv_id][i] = dna[indiv_id][j];
}


__device__
void dna_gpu_do_switch(char** dna, int indiv_id, int pos) {
    if (dna[indiv_id][pos] == '0') dna[indiv_id][pos] = '1';
    else dna[indiv_id][pos] = '0';
}

__device__
void dna_gpu_do_duplication(char** dna, char** dna_next_gen, char** dupl_seg, int indiv_id, size_t dna_size, int pos_1, int pos_2, int pos_3) {
    // Duplicate segment [pos_1; pos_2[ and insert the duplicate before pos_3
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
        seg_length = pos_2 - pos_1;
        //char seg[seg_length];
        int ix=0;

        for(int i = pos_1; i < pos_2; i++) {dupl_seg[indiv_id][ix] = dna[indiv_id][i]; ix++; }

        dna_gpu_insert(dna, dna_next_gen, indiv_id, dna_size, pos_3, dupl_seg[indiv_id],seg_length);
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
        int32_t tmp1_len = dna_size - pos_1;
        int32_t tmp2_len = pos_2;
        seg_length = tmp1_len + tmp2_len;
        //char seg[seg_length];

        int j = 0;
        for (int i = pos_1; i < dna_size; i++) {dupl_seg[indiv_id][j] = dna[indiv_id][i]; j++;}


        for (int i = 0; i < pos_2; i++) {dupl_seg[indiv_id][j] = dna[indiv_id][i]; j++;}


        dna_gpu_insert(dna, dna_next_gen, indiv_id, dna_size, pos_3, dupl_seg[indiv_id],seg_length);
    }
}

__device__
void dna_gpu_do_small_deletion(char** dna, char** dna_next_gen, int indiv_id, size_t dna_size, int pos, int nb_del) {
    // Remove promoters containing at least one nucleotide from the sequence to
    // delete

    if (pos + nb_del <= dna_size) { // the deletion does not contain the origin of
        // replication
        // Do the deletion
        dna_gpu_remove(dna, dna_next_gen, indiv_id, dna_size, pos, pos + nb_del);
    }
    else { // the deletion contains the origin of replication
        // Do the deletion
        int32_t nb_del_at_pos_0 = nb_del - dna_size + pos;

        dna_gpu_remove(dna, dna_next_gen, indiv_id, dna_size, pos, dna_size);

        dna_gpu_remove(dna, dna_next_gen, indiv_id, dna_size, 0, nb_del_at_pos_0);
    }
}

__device__
void dna_gpu_do_deletion(char** dna, char** dna_next_gen, int indiv_id, size_t dna_size, int pos_1, int pos_2) {
    // Delete segment going from pos_1 (included) to pos_2 (excluded)
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

        dna_gpu_remove(dna, dna_next_gen, indiv_id, dna_size, pos_1, pos_2);
    }
    else { // if (pos_1 >= pos_2)
        // The segment to delete includes the origin of replication.
        // The deletion process will be done in two steps.
        //
        //                                            ,->
        //    pos_2                 pos_1            |      -> 0-
        //      |                     |                   -       - pos_2
        // 0--------------------------------->     pos_1 -         -
        // =====                      =======            -         -
        //  tmp2                        tmp1              -       -
        //                                                  -----
        //
        //

        // Remove promoters containing at least one nucleotide from the sequence
        // to delete
        dna_gpu_remove(dna, dna_next_gen, indiv_id, dna_size, pos_1, dna_size); // delete tmp1 from genome
        dna_gpu_remove(dna, dna_next_gen, indiv_id, dna_size, 0, pos_2);       // delete tmp2 from genome
    }
}
