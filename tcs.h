#ifndef TCS_H_INCLUDED
#define TCS_H_INCLUDED

#include "myfunction.h"
#include "myutil.h"

#define INFINITY 1000000000
#define CORRECT_SCORE 20
#define TOO_BIG_SCORE 10
#define LWEdge_INTINT 0
#define LWEdge_HAPINT 1
#define LWEdge_HAPHAP 2

//-- input
char *inFile;
int thread_num;
int is_mask_ambiguous;
int is_merge_intermediate;
char *outFile;
int tcs_is_nucleotide_base;

//=================================

typedef struct my_iterator_t{
    AList_l *list;
    int count;
}MyIterator;

typedef struct taxaitem_t{
    int id;
    AList_l *paths;
    char *name;
    char *characters;
    int l_characters;
    int numduplicates;
    double oweight;
    double weight;
    int visited;
    int levelNumber;
    int resolved;
    int isAmbiguous;
    AList_l *realdist;
	AList_l *compdist;
    AList_l *nbor;
    AList_l *newconnections;
    AList_l *dupnames;
    AList_l *metricdist;
    void *parentComponent;
    int minRealDist;
    int isIntermediate;
    int isTip;
}TaxaItem;

typedef struct component_t{
    int id;
    AList_l *compdist;
    void *mindist;
    AList_l *taxa;
}Component;

typedef struct distance_t{
    int source;
    int destination;
    Component *sc;
    Component *dc;
    int distance;
    int marked;
}Distance;

typedef struct path_t{
    TaxaItem *source;
    TaxaItem *dest;
    int type;
    int resolved;
    AList_l *edges;
    AList_l *differences;
    int isAmbiguous;
}Path;

typedef struct lwedge_t{
    TaxaItem *source;
    TaxaItem *dest;
    int type;
    char *label;
    AList_l *third;
    int inPath;
    Path *path;
}LWEdge;

typedef struct tcs_node_t{
    int id;
    char *name;
    char *seq;
    int is_intermediate;
    int freq;
    LList_l *in_edges;
    LList_l *out_edges;
    int is_output;
    AList_l *dupnames;
}TCSNode;

typedef struct tcs_edge_t{
    TCSNode *from;
    TCSNode *to;
    int weight;
    int is_output;
}TCSEdge;



//--
int seq_num;
int seq_len;
AList_l *components;
AList_l *alltaxa;
AList_l *realtaxa;
char **maxWtaxa;
int maxParsimonyDistance;
int evaluateMetric_is_inf;
pthread_mutex_t *locker;
int *evaluateMetric_results;
AList_l *edges;

//=================================
int main_tcs(int argc, char **argv);
void load_parameters(int argc, char **argv);
void print_usage();
//=================================
void tcs_read_input_file();
void tcs_get_data(char *file, AList_l *names, AList_l *seqs);
int *tcs_get_data_mask(char *file, int l_mask);
//=================================
MyIterator *new_my_iterator(AList_l *list);
int my_iterator_has_more(MyIterator *iter);
void *my_iterator_next(MyIterator *iter);
TaxaItem *new_taxaitem1(char *n, int ident, char *chars, int l_chars);
TaxaItem *new_taxaitem2(char *n, int length, int ident);
void free_taxaitem(TaxaItem *item);
Component *new_component(int ident);
Distance *new_distance1();
Distance *new_distance2(int src, int dest, int dist);
Distance *new_distance3(int src, int dest, int dist, Component *srcc, Component *dstc);
void clone_distance(Distance *dest, Distance *src);
Path *new_path(TaxaItem *mySourcet, TaxaItem *myDestt, int myType);
LWEdge *new_LWEdge(TaxaItem *myDestt, TaxaItem *mySourcet, int myType, int label);
void free_component(Component *cpt);
//=================================
void calculate_connection_limit();
double calcPars(int j, int m, int it);
void build_matrix();
void build_matrix_one(int thread_index);
//=================================
int connectTaxa(TaxaItem *source);
int connectTaxa(TaxaItem *source);
int bestMetric(Component *sourcec, Component *destc, TaxaItem *sourcet, TaxaItem *destt, int limit, AList_l *destCandidates);
int combineComponents(Component *componentOne, Component *componentTwo);
void print_taxa(TaxaItem *item);
void print_distance(Distance *dis);
void print_component(Component *cpt);
void remove_alist_l(AList_l *list, int64_t v);
void recalcMinDistance(Component *collapse, int removeid);
void recalcDist(AList_l *allvect, TaxaItem *sourcebordert, TaxaItem *destbordert, int length);
void recalcDist_one(int64_t para);
int evaluateMetric(AList_l *allvect, int TOO_SMALL_SCORE, int numIntermediates, int limit);
void evaluateMetric_one(int64_t para);
void set_evaluateMetric_is_inf(int is_inf);
int get_evaluateMetric_is_inf();
void connectComponents(Component *sourcec, Component *destc, Distance *newdist, AList_l *candidates);
void updateDistance(Component *sourcec, Component *destc, Distance *newdist);
void addLWEdge(TaxaItem *sourcet, TaxaItem *destt);
void addIntermediate(TaxaItem *newtaxa, TaxaItem *before, TaxaItem *after, int after_distance);
int connectComponents2();
//=================================
void calculateOutgroupWeights();
void get_LWEdge_diff(LWEdge *edge, char *res_c1, char *res_c2);
void printGraph();
//=================================
void get_TCSEdge_diff(TCSEdge *edge, char *res_c1, char *res_c2);
void output_graph2();
//=================================

#endif // TCS_H_INCLUDED
