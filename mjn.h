#ifndef MJN_H_INCLUDED
#define MJN_H_INCLUDED

#include "myfunction.h"
#include "myutil.h"

#define max_vertex_num 1<<30

#define BONUS 20
#define SHORTCUTPENALTY 10
#define LONGPENALTY 5
#define TCS_MIN_SCORE -2147483648

typedef struct entry_dup_t{
    int hash;
    char *seq;
    int l_seq;
    int64_t value;
    int freq;
    AList_l *dupnames;
    int64_t next;
}Entry_dup;

typedef struct entry_edge_t{
    int64_t key1;
    int64_t key2;
    int64_t value;
    int64_t next;
}Entry_edge;

typedef struct entry_med_t{
    int64_t key1;
    int64_t key2;
    int64_t key3;
    AList_l *medes;
    AList_i *costs;
    int freq;
    int64_t next;
}Entry_med;

typedef struct sequence_t{
    char *name;
    char *seq;
    int16_t *bit_seq;
    int dup_num;
    AList_l *dupnames;
    LList_l *in_edges;
    LList_l *out_edges;
    int index;
}Sequence;

typedef struct msn_edge_t{
    Sequence *from;
    Sequence *to;
    int weight;
}Edge;

typedef struct tree_node_t{
    void *parent;
    int index;
}TNode;

typedef struct tcs_tmp_result_t{
    int index1;
    int index2;
    int compU;
    int compV;
    int dP;
    int clustDist;
    int score;
    int64_t *left;
    int64_t *right;
}TCS_tmp1;

//=============================================
char *algorithm_name;
int msn_thread_num;
int epsilon;
char *msn_inFile;
char *msn_outFile;
ThreadPool *msn_thread_pool;
pthread_mutex_t *lock_for_multi_thread;
pthread_cond_t *cond_for_multi_thread;
int is_nucleotide_base;
//--
AList_i *tcs_component_ids;
int *path_lengths;
int64_t r_path_lengths;
//--
int n_seq;
int l_seq;
int n_seq_bit;
int l_seq_bit;
AList_l *data;
int *weigths;
int *distances;
int *pre_distances;
int r_distance;
int pre_r_distance;
int min_distance;
int max_distance;
int *int_2_distance;
int *base_2_int;
int is_use_avx;
int is_use_bit_operation_for_distance;
//--
int l_table_verteres;
int64_t *table_verteres;
LList_l *all_verteres;
int l_table_edges;
int64_t *table_edges;
LList_l *all_edges;
//--
AList_l *all_meds;
int l_table_medes;
int64_t *table_medes;
//--
int msn_intermediate_index;
int mjn_compute_number;
int msn_compute_number;



//=============================================
int main_mjn(int argc, char **argv);
void load_parameters2(int argc, char **argv);
void print_usage2();
//=============================================
int32_t hash_mem(char *key, int l_key);
Entry_dup *put_hash_seqs(int64_t *table, int l_table, char *seq, int l_seq, int64_t value);
Entry_dup *get_hash_seqs(int64_t *table, int l_table, char *seq, int l_seq);
void remove_hash_seqs(int64_t *table, int l_table, int64_t value);
Entry_edge *put_hash_edges(int64_t *table, int l_table, int64_t key1, int64_t key2, int64_t value);
Entry_edge *get_hash_edges(int64_t *table, int l_table, int64_t key1, int64_t key2);
void remove_hash_edges(int64_t *table, int l_table, int64_t key1, int64_t key2);
Entry_med *put_hash_meds(int64_t *table, int l_table, int64_t key1, int64_t key2, int64_t key3, AList_l *medes, AList_i *costs);
void sort_triple(int64_t key1, int64_t key2, int64_t key3, int64_t *res1, int64_t *res2, int64_t *res3);
//=============================================
Sequence *add_vertex0(char *name, char *seq, int dup_num, AList_l *dupnames);
Sequence *add_vertex(char *name, char *seq, int dup_num);
Edge *add_edge(Sequence *from, Sequence *to, int weight);
Edge *get_edge(Sequence *from, Sequence *to);
AList_l *get_all_edges(Sequence *v);
Sequence *get_edge_opposite(Edge *edge, Sequence *v);
void remove_vertex(Sequence *v);
void clear_all_edges();
void clear_all_edges_one(int64_t para);
void free_entry_edge(Entry_edge *e);
//=============================================
void msn_read_data();
int *get_data_mask(int l_mask);
void is_site_i_j_is_same(int64_t para);
void calculate_distances();
void calculate_distances_one(int64_t para);
int get_pairwise_distance(Sequence *seq1, Sequence *seq2, int16_t *buff);
void set_distance(int i, int j, int dis);
int get_distance(int i, int j);
void calculate_int_2_distance();
int16_t *convert_base_seq_2_bit_seq(char *seq);
//=============================================
void compute_modified_tcs();
int findIntermediates(Pair *intPair, Sequence *u, Sequence *v, int dist);
void newCompositePath(Sequence *start, Sequence *end, int dist);
void newCompositePath_one(int64_t para);
void computeScore_one(int64_t para);
int computeScore(int index1, int index2, int compU, int compV, int dP, int clustDist);
int get_path_length(int index1, int index2);
void set_path_length(int index1, int index2, int dis);
void ensure_path_lengths();
void update_path_lengths();
void update_path_lengths_with_edge(int u_index, int v_index, int weight);
void update_path_lengths_with_edge_one1(int64_t para);
void update_path_lengths_with_edge_one2(int64_t para);
//=============================================
void compute_msn_only();
//=============================================
void compute_mjn_graph();
void compute_msn(AList_l *feasibleLinks);
void compute_mjn();
void calculate_min_cost_one(int64_t para);
void calculate_min_cost_one2(int64_t para);
void computeQuasiMedianSeqs(char *seqA, char *seqB, char *seqC, SBuilder *sb, int len, AList_l *list, AList_i *costs);
int computeCost(char *seqA, char *seqB, char *seqC, char *med);
int pairwiseDistance(char *s1, char *s2);
int remove_obsolete_vertexes();
//=============================================
void msn_output_graph();
//=============================================

#endif // MJN_H_INCLUDED
