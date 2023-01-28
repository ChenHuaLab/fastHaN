#include "mjn.h"

//=============================================

int main_mjn(int argc, char **argv){
    load_parameters2(argc, argv);

    mylog("starting...");

    int i, j;

    msn_thread_pool=new_thread_pool(msn_thread_num);
    lock_for_multi_thread=my_new(1, sizeof(pthread_mutex_t));
    cond_for_multi_thread=my_new(1, sizeof(pthread_cond_t));
    pthread_mutex_init(lock_for_multi_thread, NULL);
    pthread_cond_init(cond_for_multi_thread, NULL);

    msn_read_data();
    if(is_use_bit_operation_for_distance) calculate_int_2_distance();

    l_table_verteres=n_seq*100;
    table_verteres=my_new(l_table_verteres, sizeof(int64_t));
    all_verteres=new_llist_l();
    l_table_edges=n_seq*n_seq;
    table_edges=my_new(l_table_edges, sizeof(int64_t));
    all_edges=new_llist_l();

    for(i=0;i<n_seq;i++){
        Sequence *seq=data->elementData[i];
        add_vertex0(seq->name, seq->seq, seq->dup_num, seq->dupnames);
    }

    if(strcmp(algorithm_name, "modified_tcs")==0) compute_modified_tcs();
    else if(strcmp(algorithm_name, "msn")==0) compute_msn_only();
    else compute_mjn_graph();

    mylog("----------");
    mylog("writing graph into file...");
    msn_output_graph();

    free_thread_pool(msn_thread_pool);
    mylog("finished!");

    return 0;
}

void load_parameters2(int argc, char **argv){
    algorithm_name=argv[0];
    argc--;
    argv++;

    if(argc<4 || argc>8) print_usage2();

    msn_thread_num=8;
    epsilon=0;
    msn_intermediate_index=1;
    mjn_compute_number=1;
    msn_compute_number=1;
    distances=NULL;

    is_use_avx=1;
    is_use_bit_operation_for_distance=1;

    r_path_lengths=0;
    path_lengths=NULL;

    int i;
    for(i=0;i<argc;i++){
        if(strcmp(argv[i], "-i")==0){
            i++;
            msn_inFile=argv[i];
        }else if(strcmp(argv[i], "-t")==0){
            i++;
            msn_thread_num=atoi(argv[i]);
        }else if(strcmp(argv[i], "-e")==0){
            i++;
            epsilon=atoi(argv[i]);
        }else if(strcmp(argv[i], "-o")==0){
            i++;
            msn_outFile=argv[i];
        }else{
            fprintf(stderr, "'%s' parameter is not exist, please check!\n");
            exit(0);
        }
    }

    if(msn_thread_num<1) msn_thread_num=1;
    if(msn_thread_num>32767) msn_thread_num=32767;
}

void print_usage2(){
    printf("===========================\n");
    printf("Author:   ChiLianjiang\n");
    printf("E-mail:   chilianjiang@126.com\n");
    printf("Date:     2021-04-13\n");
    printf("Version:  1.0\n\n");
    printf("Usage: fastHapNetwork %s [arguments]\n", algorithm_name);
    printf("  input:\n");
    printf("    -i    input phylip format file\n");
    printf("  options:\n");
    if(strcmp(algorithm_name, "modified_tcs")==0 || strcmp(algorithm_name, "mjn")==0)
    printf("    -t    (int)thread number(default:8)\n");
    if(strcmp(algorithm_name, "msn")==0 || strcmp(algorithm_name, "mjn")==0)
    printf("    -e    (int)epsilon(default:0)\n");
    printf("  output:\n");
    printf("    -o    output prefix(graph file)(GML and json format)\n\n");
    exit(0);
}

//=============================================

int32_t hash_mem(char *key, int l_key){
    int32_t i, h=0;
    for(i=0;i<l_key;i++) h=31*h+key[i];
    return h;
}

Entry_dup *put_hash_seqs(int64_t *table, int l_table, char *seq, int l_seq, int64_t value){
    int hash=hash_mem(seq, l_seq);
    int h_index=((unsigned int)hash)%l_table;
    if(h_index<0){fprintf(stderr, "error hash!\n");exit(0);}

    Entry_dup *e=table[h_index];
    while(e){
        if(e->hash==hash && e->l_seq==l_seq && memcmp(e->seq, seq, l_seq)==0){
            e->freq++;
            return e;
        }
        e=e->next;
    }
    e=my_new(1, sizeof(Entry_dup));
    e->hash=hash;
    e->seq=str_copy_with_len(seq, l_seq);
    e->l_seq=l_seq;
    e->value=value;
    e->freq=1;
    e->next=table[h_index];
    table[h_index]=e;
    return e;
}

Entry_dup *get_hash_seqs(int64_t *table, int l_table, char *seq, int l_seq){
    int hash=hash_mem(seq, l_seq);
    int h_index=((unsigned int)hash)%l_table;
    if(h_index<0){fprintf(stderr, "error hash!\n");exit(0);}

    Entry_dup *e=table[h_index];
    while(e){
        if(e->hash==hash && e->l_seq==l_seq && memcmp(e->seq, seq, l_seq)==0){
            e->freq++;
            return e;
        }
        e=e->next;
    }
    return NULL;
}

void remove_hash_seqs(int64_t *table, int l_table, int64_t value){
    Sequence *seq=(Sequence *)value;

    int hash=hash_mem(seq->seq, l_seq);
    int h_index=((unsigned int)hash)%l_table;
    if(h_index<0){fprintf(stderr, "error hash!\n");exit(0);}

    Entry_dup *e=table[h_index];
    Entry_dup *pre=NULL;

    for(;e;e=e->next){
        if(e->value==value){
            if(pre==NULL) table[h_index]=e->next;
            else pre->next=e->next;
            //--
            free(e->seq);
            free(e);
            //--
            return;
        }
        pre=e;
    }
}

Entry_edge *put_hash_edges(int64_t *table, int l_table, int64_t key1, int64_t key2, int64_t value){
    int h_index=(int)(((unsigned)(key1^key2))%((int64_t)l_table));
    if(h_index<0){fprintf(stderr, "error hash!\n");exit(0);}

    Entry_edge *e=table[h_index];
    while(e){
        if(e->key1==key1 && e->key2==key2) return e;
        e=e->next;
    }
    e=my_new(1, sizeof(Entry_edge));
    e->key1=key1;
    e->key2=key2;
    e->value=value;
    e->next=table[h_index];
    table[h_index]=e;
    return e;
}

Entry_edge *get_hash_edges(int64_t *table, int l_table, int64_t key1, int64_t key2){
    int h_index=(int)(((unsigned)(key1^key2))%((int64_t)l_table));
    if(h_index<0){fprintf(stderr, "error hash!\n");exit(0);}

    Entry_edge *e=table[h_index];
    while(e){
        if(e->key1==key1 && e->key2==key2) return e;
        e=e->next;
    }
    return NULL;
}

void remove_hash_edges(int64_t *table, int l_table, int64_t key1, int64_t key2){
    int h_index=(int)(((unsigned)(key1^key2))%((int64_t)l_table));
    if(h_index<0){fprintf(stderr, "error hash!\n");exit(0);}

    Entry_edge *e=table[h_index];
    Entry_edge *pre=NULL;
    for(;e;e=e->next){
        if(e->key1==key1 && e->key2==key2){
            if(pre==NULL) table[h_index]=e->next;
            else pre->next=e->next;
            free(e);
            return;
        }
        pre=e;
    }
}

Entry_med *put_hash_meds(int64_t *table, int l_table, int64_t key1, int64_t key2, int64_t key3, AList_l *medes, AList_i *costs){
    //sort_triple(key1, key2, key3, &key1, &key2, &key3);
    int h_index=(int)(((unsigned)(key1^key2^key3))%((int64_t)l_table));
    if(h_index<0){fprintf(stderr, "error hash!\n");exit(0);}

    Entry_med *e=table[h_index];
    while(e){
        if(e->key1==key1 && e->key2==key2 && e->key3==key3){e->freq++;return e;}
        e=e->next;
    }
    e=my_new(1, sizeof(Entry_med));
    e->key1=key1;
    e->key2=key2;
    e->key3=key3;
    e->medes=medes;
    e->costs=costs;
    e->freq=1;
    e->next=table[h_index];
    table[h_index]=e;
    return e;
}

void sort_triple(int64_t key1, int64_t key2, int64_t key3, int64_t *res1, int64_t *res2, int64_t *res3){
    if(key1>=key2 && key1>=key3 && key2>=key3){*(res1)=key1;*(res2)=key2;*(res3)=key3;}
    else if(key1>=key2 && key1>=key3 && key3>=key2){*(res1)=key1;*(res2)=key3;*(res3)=key2;}
    else if(key2>=key1 && key2>=key3 && key1>=key3){*(res1)=key2;*(res2)=key1;*(res3)=key3;}
    else if(key2>=key1 && key2>=key3 && key3>=key1){*(res1)=key2;*(res2)=key3;*(res3)=key1;}
    else if(key3>=key1 && key3>=key2 && key1>=key2){*(res1)=key3;*(res2)=key1;*(res3)=key2;}
    else {*(res1)=key3;*(res2)=key2;*(res3)=key1;}
}

//=============================================

Sequence *add_vertex0(char *name, char *seq, int dup_num, AList_l *dupnames){
    Sequence *v=my_new(1, sizeof(Sequence));
    Entry_dup *e=put_hash_seqs(table_verteres, l_table_verteres, seq, l_seq, v);
    if(e->freq!=1){
        fprintf(stderr, "[error]already exists seq=%s\n", seq);
        exit(0);
    }
    v->name=str_copy(name);
    v->seq=str_copy_with_len(seq, l_seq);
    if(is_use_bit_operation_for_distance) v->bit_seq=convert_base_seq_2_bit_seq(seq);
    v->dup_num=dup_num;
    v->in_edges=new_llist_l();
    v->out_edges=new_llist_l();
    v->dupnames=dupnames;
    llist_l_add(all_verteres, v);
    if(all_verteres->size>max_vertex_num){fprintf(stderr, "cannot process so much vertexes: %d\n", all_verteres->size);exit(0);}
    return v;
}

Sequence *add_vertex(char *name, char *seq, int dup_num){
    Sequence *v=my_new(1, sizeof(Sequence));
    Entry_dup *e=put_hash_seqs(table_verteres, l_table_verteres, seq, l_seq, v);
    if(e->freq!=1){
        fprintf(stderr, "[error]already exists seq=%s\n", seq);
        exit(0);
    }
    v->name=str_copy(name);
    v->seq=str_copy_with_len(seq, l_seq);
    if(is_use_bit_operation_for_distance) v->bit_seq=convert_base_seq_2_bit_seq(seq);
    v->dup_num=dup_num;
    v->in_edges=new_llist_l();
    v->out_edges=new_llist_l();
    v->dupnames=new_alist_l(1);
    alist_l_add(v->dupnames, str_copy(name));
    llist_l_add(all_verteres, v);
    if(all_verteres->size>max_vertex_num){fprintf(stderr, "cannot process so much vertexes: %d\n", all_verteres->size);exit(0);}
    return v;
}

Edge *add_edge(Sequence *from, Sequence *to, int weight){
    if(from==to){
        fprintf(stderr, "[error]connect edge to itself , %s-%s\n", from->name, to->name);
        exit(0);
    }

    Edge *edge=get_edge(from, to);
    if(edge) return edge;

    edge=get_edge(to, from);
    if(edge){
        fprintf(stderr, "[warnning]connect circle edge, %s-%s\n", from->name, to->name);
        return edge;
    }

    edge=my_new(1, sizeof(Edge));
    edge->from=from;
    edge->to=to;
    edge->weight=weight;
    Entry_edge *entry=put_hash_edges(table_edges, l_table_edges, from, to, edge);
    llist_l_add(from->out_edges, edge);
    llist_l_add(to->in_edges, edge);
    llist_l_add(all_edges, edge);

    return edge;
}

inline Edge *get_edge(Sequence *from, Sequence *to){
    Entry_edge *entry=get_hash_edges(table_edges, l_table_edges, from, to);
    return entry ? entry->value:NULL;
}

AList_l *get_all_edges(Sequence *v){
    AList_l *list=new_alist_l(v->in_edges->size+v->out_edges->size);
    LListNode_l *n;
    for(n=v->in_edges->first;n;n=n->next) alist_l_add(list, n->value);
    for(n=v->out_edges->first;n;n=n->next) alist_l_add(list, n->value);
    return list;
}

inline Sequence *get_edge_opposite(Edge *edge, Sequence *v){
    return edge->from==v ? (edge->to):(edge->from);
}

void remove_vertex(Sequence *seq){
    LListNode_l *n;
    for(n=seq->in_edges->first;n;n=n->next){
        Edge *e=n->value;
        llist_l_remove2(e->from->out_edges, e);
        llist_l_remove2(all_edges, e);
        remove_hash_edges(table_edges, l_table_edges, e->from, e->to);
        free(e);
    }
    for(n=seq->out_edges->first;n;n=n->next){
        Edge *e=n->value;
        llist_l_remove2(e->to->in_edges, e);
        llist_l_remove2(all_edges, e);
        remove_hash_edges(table_edges, l_table_edges, e->from, e->to);
        free(e);
    }
    //--
    llist_l_remove2(all_verteres, seq);
    remove_hash_seqs(table_verteres, l_table_verteres, seq);
    free(seq->name);
    free(seq->seq);
    if(is_use_bit_operation_for_distance) free_align_ptr(seq->bit_seq);
    free_llist_l(seq->in_edges);
    free_llist_l(seq->out_edges);
    free(seq);
}

void clear_all_edges(){
    LListNode_l *n=all_edges->first;
    while(n){
        LListNode_l *c=n;
        n=(LListNode_l *)n->next;
        free(c->value);
        free(c);
    }
    free(all_edges);
    all_edges=new_llist_l();
    //--
    split_span_for_threads(msn_thread_num, 0, l_table_edges);
    int i;
    for(i=0;i<msn_thread_num;i++) thread_pool_add_worker(msn_thread_pool, clear_all_edges_one, i);
    thread_pool_invoke_all(msn_thread_pool);
    memset(table_edges, 0, l_table_edges*sizeof(int64_t));
    //--
    for(n=all_verteres->first;n;n=n->next){
        Sequence *v=n->value;
        free_llist_l(v->in_edges);
        v->in_edges=new_llist_l();
        free_llist_l(v->out_edges);
        v->out_edges=new_llist_l();
    }
}

void clear_all_edges_one(int64_t para){
    int i=(int)(para&0xFFFFFFFFL);
    int start=starts_thread[i];
    int end=ends_thread[i];

    for(i=start;i<end;i++) free_entry_edge(table_edges[i]);
}

void free_entry_edge(Entry_edge *e){
    if(!e) return;
    free_entry_edge(e->next);
    free(e);
}

//=============================================

Hash_li *global_samePosAs;
char *global_all_seq;

void msn_read_data(){
    sprintf(loginfo, "reading data file(%s)...", msn_inFile);mylog(loginfo);

    int i, j, k, t, line_len, max_len=1000000000;

    data=new_alist_l(16);

    GzStream *gz1=gz_stream_open(msn_inFile, "r");

    //-- first line
    gz_read_util(gz1, '\n', line, max_len, &line_len);
    line_len=chmop_with_len(line, line_len);
    str_replace_char_with_no_copy(line, ' ', '\t');
    i=0;
    while(line[++i]!='\t');
    line[i]='\0';
    n_seq=atoi(line);
    l_seq=atoi(line+i+1);

    int *mask=get_data_mask(l_seq);
    int new_l_seq=0;

    int l_table=2*n_seq;
    int64_t *table=my_new(l_table, sizeof(int64_t));

    while(gz_read_util(gz1, '\n', line, max_len, &line_len)){
        line_len=chmop_with_len(line, line_len);
        str_replace_char_with_no_copy(line, '\t', ' ');
        //--
        i=0;
        while(line[++i]!=' ');
        line[i]='\0';
        char *name=line;
        int l_name=i;
        //--
        while(line[++i]==' ');
        char *seq=line+i;
        //--
        SBuilder *sb=new_s_builder(l_seq);
        for(i=0;i<l_seq;i++){
            if(!mask[i]) s_builder_add_char(sb, toupper(seq[i]));
        }
        //--
        Entry_dup *e=put_hash_seqs(table, l_table, sb->str, sb->size, 0);
        if(e->freq==1){
            e->dupnames=new_alist_l(2);
            alist_l_add(e->dupnames, str_copy_with_len(name, l_name));
            alist_l_add(data, e);
        }else alist_l_add(e->dupnames, str_copy_with_len(name, l_name));
        //--
        new_l_seq=sb->size;
        free_s_builder(sb);
    }
    gz_stream_destory(gz1);
    free(mask);
    l_seq=new_l_seq;

    sprintf(loginfo, "seq_num=%d\tmark_dup_num=%d", n_seq, data->size);mylog(loginfo);
    n_seq=data->size;

    for(i=0;i<data->size;i++){
        Entry_dup *e=data->elementData[i];
        Sequence *s=my_new(1, sizeof(Sequence));
        s->name=e->dupnames->elementData[0];
        s->seq=e->seq;
        s->dup_num=e->freq;
        s->dupnames=e->dupnames;
        data->elementData[i]=s;
        free(e);
    }
    free(table);

    int *samePosAs=my_new(l_seq, sizeof(int));
    for(i=0;i<l_seq;i++) samePosAs[i]=i;

    //-- transpose data matrix
    global_all_seq=my_new((int64_t)n_seq*(int64_t)l_seq, sizeof(char));
    for(i=0;i<n_seq;i++){
        char *seq=((Entry_dup *)(data->elementData[i]))->seq;
        for(j=0;j<l_seq;j++) global_all_seq[(int64_t)j*(int64_t)n_seq+(int64_t)i]=seq[j];
    }

    if(msn_thread_num<2){
        char *i2j=my_new(256, sizeof(char));
        char *j2i=my_new(256, sizeof(char));
        for(i=0;i<l_seq;i++){
            char *seq1=global_all_seq+(int64_t)i*(int64_t)n_seq;
            for(j=i+1;j<l_seq;j++){
                memset(i2j, 0, 256*sizeof(char));
                memset(j2i, 0, 256*sizeof(char));
                char *seq2=global_all_seq+(int64_t)j*(int64_t)n_seq;
                int same=1;
                for(t=0;same&&t<n_seq;t++){
                    char chari=seq1[t];
                    char charj=seq2[t];
                    //--
                    if(i2j[chari]==(char)0){
                        i2j[chari]=charj;
                        if (j2i[charj]==(char)0) j2i[charj]=chari;
                        else if (j2i[charj]!=chari) same=0;
                    }else if(i2j[chari]!=charj) same=0;
                }
                if(same){
                    samePosAs[j]=samePosAs[i];
                    break;
                }
            }
        }
        free(i2j);
        free(j2i);
    }else{
        global_samePosAs=new_hash_li1(l_seq);
        for(i=0;i<l_seq;i++) thread_pool_add_worker(msn_thread_pool, is_site_i_j_is_same, i);
        thread_pool_invoke_all(msn_thread_pool);

        for(i=0;i<l_seq;i++){
            for(j=i+1;j<l_seq;j++){
                int64_t key=(int64_t)i*(int64_t)l_seq+(int64_t)+j;
                if(hash_li_get(global_samePosAs, key)){
                    samePosAs[j]=samePosAs[i];
                    break;
                }
            }
        }
        free_hash_li(global_samePosAs);
    }

    SBuilder **all_bases=my_new(n_seq, sizeof(SBuilder *));
    for(i=0;i<n_seq;i++) all_bases[i]=new_s_builder(16);

    new_l_seq=0;
    int *flag=my_new(l_seq, sizeof(int));

    for(i=0;i<l_seq;i++){
        char *seq=global_all_seq+(int64_t)i*(int64_t)n_seq;
        if(samePosAs[i]==i){
            new_l_seq++;
            for(j=0;j<n_seq;j++) s_builder_add_char(all_bases[j], seq[j]);
        }
        flag[samePosAs[i]]++;
    }
    free(samePosAs);
    free(global_all_seq);

    weigths=my_new(new_l_seq, sizeof(int));
    for(i=0,j=0;i<l_seq;i++){
        if(flag[i]) weigths[j++]=flag[i];
    }
    free(flag);

    sprintf(loginfo, "seq_len=%d\tmark_dup_len=%d", l_seq, new_l_seq);mylog(loginfo);
    l_seq=new_l_seq;

    for(i=0;i<n_seq;i++){
        Sequence *seq=data->elementData[i];
        free(seq->seq);
        seq->seq=all_bases[i]->str;
        free(all_bases[i]);
    }
    free(all_bases);

    if(is_use_bit_operation_for_distance){
        if(is_use_avx){
            n_seq_bit=l_seq/128;
            if((n_seq_bit*128)<l_seq) n_seq_bit++;
            l_seq_bit=n_seq_bit*(256/16);
            sprintf(loginfo, "l_seq=%d\tnum_seq_bit=%d\tlen_seq_bit=%d", l_seq, n_seq_bit, l_seq_bit);mylog(loginfo);
        }else{
            l_seq_bit=l_seq/8;
            if((l_seq_bit*8)<l_seq) l_seq_bit++;
            sprintf(loginfo, "l_seq=%d\tl_seq_bit=%d", l_seq, l_seq_bit);mylog(loginfo);
        }
    }
}

int *get_data_mask(int l_mask){
    int *mask=my_new(l_mask, sizeof(int));

    int i, line_len, max_len=1000000000;

    is_nucleotide_base=1;

    GzStream *gz1=gz_stream_open(msn_inFile, "r");
    //-- first line
    gz_read_util(gz1, '\n', line, max_len, &line_len);
    while(gz_read_util(gz1, '\n', line, max_len, &line_len)){
        line_len=chmop_with_len(line, line_len);
        str_replace_char_with_no_copy(line, '\t', ' ');
        //--
        i=0;
        while(line[++i]!=' ');line[i]='\0';
        while(line[++i]==' ');
        char *seq=line+i;
        //--
        for(i=0;i<l_seq;i++){
            char c=toupper(seq[i]);
            if (c == '-' || c == '.' ||
                c == 'R' || c == 'M' || c == 'W' || c == 'S' || c == 'K' || c == 'Y' ||
                c == 'H' || c == 'V' || c == 'D' || c == 'B' || c == 'X' || c == 'N'){
                c='?';
            }
            //--
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != '?'){
                is_nucleotide_base=0;
                goto pos;
            }
        }
    }
    pos:
    gz_stream_destory(gz1);

    gz1=gz_stream_open(msn_inFile, "r");
    //-- first line
    gz_read_util(gz1, '\n', line, max_len, &line_len);
    while(gz_read_util(gz1, '\n', line, max_len, &line_len)){
        line_len=chmop_with_len(line, line_len);
        str_replace_char_with_no_copy(line, '\t', ' ');
        //--
        i=0;
        while(line[++i]!=' ');line[i]='\0';
        while(line[++i]==' ');
        char *seq=line+i;
        //--
        for(i=0;i<l_seq;i++){
            char c=toupper(seq[i]);
            if(is_nucleotide_base){
                if (c == '-' || c == '.' ||
                    c == 'R' || c == 'M' || c == 'W' || c == 'S' || c == 'K' || c == 'Y' ||
                    c == 'H' || c == 'V' || c == 'D' || c == 'B' || c == 'X' || c == 'N'){
                    mask[i]=1;
                }
            }else{
                if (c == '-' || c == '.') mask[i]=1;
            }
        }
    }
    gz_stream_destory(gz1);

    int ambiguous_sites_num=0;
    for(i=0;i<l_mask;i++){
        if(mask[i]) ambiguous_sites_num++;
    }
    sprintf(loginfo, "ambiguous_sites_num=%d", ambiguous_sites_num);mylog(loginfo);

    return mask;
}

void is_site_i_j_is_same(int64_t para){
    int j, i=(int)(para&0xFFFFFFFFL);

    char *i2j=my_new(256, sizeof(char));
    char *j2i=my_new(256, sizeof(char));

    char *seq1=global_all_seq+(int64_t)i*(int64_t)n_seq;
    for(j=i+1;j<l_seq;j++){
        memset(i2j, 0, 256*sizeof(char));
        memset(j2i, 0, 256*sizeof(char));
        char *seq2=global_all_seq+(int64_t)j*(int64_t)n_seq;

        int t, same=1;
        for(t=0;same&&t<n_seq;t++){
            char chari=seq1[t];
            char charj=seq2[t];
            //--
            if(i2j[chari]==(char)0){
                i2j[chari]=charj;
                if (j2i[charj]==(char)0) j2i[charj]=chari;
                else if (j2i[charj]!=chari) same=0;
            }else if(i2j[chari]!=charj) same=0;
        }
        if(same){
            int64_t key=(int64_t)i*(int64_t)l_seq+(int64_t)j;
            pthread_mutex_lock(lock_for_multi_thread);
            hash_li_put(global_samePosAs, key, 1);
            pthread_mutex_unlock(lock_for_multi_thread);
            break;
        }
    }

    free(i2j);
    free(j2i);
}

void calculate_distances(){
    mylog("calculating distances...");

    pre_distances=distances;
    pre_r_distance=r_distance;

    r_distance=all_verteres->size;
    distances=my_new((int64_t)r_distance*(int64_t)r_distance, sizeof(int));

    min_distance=1<<30;
    max_distance=-1;

    int i;
    split_span_for_threads(msn_thread_num, 0, all_verteres->size);
    for(i=0;i<msn_thread_num;i++) thread_pool_add_worker(msn_thread_pool, calculate_distances_one, i);
    thread_pool_invoke_all(msn_thread_pool);

    if(pre_distances) free(pre_distances);

    sprintf(loginfo, "min_distance=%d\tmax_distance=%d", min_distance, max_distance);mylog(loginfo);
}

void calculate_distances_one(int64_t para){
    int thread_index=(int)(para&0xFFFFFFFFL);
    int start=starts_thread[thread_index];
    int end=ends_thread[thread_index];

    int tmp_min_distance=1<<30;
    int tmp_max_distance=-1;

    int i, j, k, num;
    LListNode_l *n1;
    LListNode_l *n2;

    int16_t *buff=NULL;
    if(is_use_avx) buff=calloc_align_ptr(l_seq_bit*sizeof(int16_t), 64);

    for(i=0,n1=all_verteres->first;n1;i++,n1=n1->next){
        if(i<start) continue;
        if(i>=end) break;
        //--
        Sequence *seq1=n1->value;
        seq1->index=i;
        for(j=0,n2=all_verteres->first;n2&&j<i;j++,n2=n2->next){
            Sequence *seq2=n2->value;
            if(i<n_seq && j<n_seq && pre_distances){
                num=pre_distances[(int64_t)i*(int64_t)pre_r_distance+(int64_t)j];
            }else{
                num=get_pairwise_distance(seq1, seq2, buff);
            }
            distances[(int64_t)i*(int64_t)r_distance+(int64_t)j]=num;
            distances[(int64_t)j*(int64_t)r_distance+(int64_t)i]=num;
            if(tmp_min_distance>num) tmp_min_distance=num;
            if(tmp_max_distance<num) tmp_max_distance=num;
        }
    }

    if(is_use_avx) free_align_ptr(buff);

    pthread_mutex_lock(lock_for_multi_thread);
    if(min_distance>tmp_min_distance) min_distance=tmp_min_distance;
    if(max_distance<tmp_max_distance) max_distance=tmp_max_distance;
    pthread_mutex_unlock(lock_for_multi_thread);
}

int get_pairwise_distance(Sequence *seq1, Sequence *seq2, int16_t *buff){
    int i, j, num=0;

    if(is_use_bit_operation_for_distance){
        int16_t *s1=seq1->bit_seq;
        int16_t *s2=seq2->bit_seq;
        if(is_use_avx){
            uint16_t *copy=buff;
            for(i=0;i<n_seq_bit;i++){
                _mm256_store_si256(buff, _mm256_xor_si256(_mm256_load_si256(s1), _mm256_load_si256(s2)));
                s1+=16;
                s2+=16;
                buff+=16;
            }
            for(i=0;i<l_seq_bit;i++) num+=int_2_distance[(int64_t)copy[i]*(int64_t)l_seq_bit+(int64_t)i];
        }else{
            for(i=0;i<l_seq_bit;i++){
                int a=((uint16_t)s1[i])^((uint16_t)s2[i]);
                num+=int_2_distance[(int64_t)a*(int64_t)l_seq_bit+(int64_t)i];
            }
        }
    }else{
        char *s1=seq1->seq;
        char *s2=seq2->seq;
        for(i=0;i<l_seq;i++){
            char c1=s1[i];
            char c2=s2[i];
            if(c1=='?' || c2=='?') continue;
            if(c1!=c2) num+=weigths[i];
        }
    }
    return num;
}

inline void set_distance(int i, int j, int dis){
    if(i<j) distances[(int64_t)i*(int64_t)r_distance+(int64_t)j]=dis;
    else distances[(int64_t)j*(int64_t)r_distance+(int64_t)i]=dis;
}

inline int get_distance(int i, int j){
    return i==j ? 0:(i<j ? (distances[(int64_t)i*(int64_t)r_distance+(int64_t)j]):(distances[(int64_t)j*(int64_t)r_distance+(int64_t)i]));
}

void calculate_int_2_distance(){
    mylog("calculating int to distance...");
    int i, j, k, p, num, max=1<<16;

    base_2_int=my_new(256, sizeof(int));
    for(i=0;i<256;i++) base_2_int[i]=4;
    base_2_int['A']=0;
    base_2_int['a']=0;
    base_2_int['C']=1;
    base_2_int['c']=1;
    base_2_int['G']=2;
    base_2_int['g']=2;
    base_2_int['T']=3;
    base_2_int['t']=3;

    int_2_distance=my_new((int64_t)max*(int64_t)l_seq_bit, sizeof(int));

    for(i=0;i<max;i++){
        for(j=0;j<l_seq_bit;j++){
            num=0;
            for(k=0;k<8;k++){
                p=j*8+k;
                if(p>=l_seq) break;
                if((i>>(k*2))&3) num+=weigths[p];
            }
            int_2_distance[(int64_t)i*(int64_t)l_seq_bit+(int64_t)j]=num;
        }
    }
}

int16_t *convert_base_seq_2_bit_seq(char *seq){
    int16_t *ptr=calloc_align_ptr(l_seq_bit*sizeof(int16_t), 32);
    int16_t *copy=ptr;

    int i, j;
    ptr--;
    for(i=0;i<l_seq;i++){
        j=i%8;
        if(j==0) ptr++;
        *(ptr)|=(base_2_int[seq[i]]<<(j*2));
    }

    return copy;
}

//=============================================

void compute_modified_tcs(){
    mylog("calculating tcs...");

    int i, j, k;
    LListNode_l *n, *n1, *n2;

    calculate_distances();

    int seqCount=all_verteres->size;
    int ncomps=all_verteres->size;

    tcs_component_ids=new_alist_i(seqCount);
    int64_t *dist2pairs=my_new(max_distance+1, sizeof(int64_t));

    for(i=0,n1=all_verteres->first;n1;i++,n1=n1->next){
        alist_i_add(tcs_component_ids, i);
        Sequence *v1=n1->value;
        for(j=0,n2=all_verteres->first;n2&&j<i;j++,n2=n2->next){
            Sequence *v2=n2->value;
            int dis=distances[(int64_t)i*(int64_t)r_distance+(int64_t)j];
            AList_l *list=dist2pairs[dis];
            if(list==NULL){
                list=new_alist_l(16);
                dist2pairs[dis]=list;
            }
            Pair *p=my_new(1, sizeof(Pair));
            p->left=v1;
            p->right=v2;
            alist_l_add(list, p);
        }
    }

    update_path_lengths();

    int M;
    for(M=1;M<=max_distance;M++){
        AList_l *vcptr=dist2pairs[M];
        if(vcptr==NULL) continue;

        sprintf(loginfo, "dis=%d\tsize=%d", M, vcptr->size);mylog(loginfo);

        int compA=-1, compB=-1;
        AList_l *other_pairs=new_alist_l(16);

        for(i=0;i<vcptr->size;i++){
            Pair *pair=vcptr->elementData[i];
            //--
            Sequence *u=pair->left;
            Sequence *v=pair->right;
            int compU=tcs_component_ids->elementData[u->index];
            int compV=tcs_component_ids->elementData[v->index];
            if(compU==compV) goto next1;
            if(compU>compV){
                int tmp1=compU;
                compU=compV;
                compV=tmp1;
                Sequence *tmp2=u;
                u=v;
                v=tmp2;
            }
            if(compA<0){
                compA=compU;
                compB=compV;
            }
            if(compU==compA && compV==compB){
                if(M==1){
                    add_edge(u, v, 1);
                    update_path_lengths_with_edge(u->index, v->index, 1);
                }else{
                    Pair *intermediates=my_new(1, sizeof(Pair));
                    int newPathLength=findIntermediates(intermediates, u, v, M);
                    int existingPath=get_path_length(((Sequence *)(intermediates->left))->index, ((Sequence *)(intermediates->right))->index);
                    int pathExists=existingPath==-1 ? 0:1;
                    if(pathExists && existingPath<newPathLength){
                        fprintf(stderr, "existingPath=%d\tnewPathLength=%d\tM=%d\n", existingPath, newPathLength, M);
                        fprintf(stderr, "[error]Shorter path already exists between these vertices!\n");
                        exit(0);
                    }else if(!pathExists || existingPath>newPathLength){
                        newCompositePath(intermediates->left, intermediates->right, newPathLength);
                    }
                    free(intermediates);
                }
            }else{
                pair->left=u;
                pair->right=v;
                alist_l_add(other_pairs, pair);
                goto next2;
            }
            //--
            next1:
            free(pair);
            next2:
            continue;
        }

        if(compA>=0){
            for(j=0;j<tcs_component_ids->size;j++){
                if(tcs_component_ids->elementData[j]< 0 || tcs_component_ids->elementData[j]==compB) tcs_component_ids->elementData[j]=compA;
                else if (tcs_component_ids->elementData[j]>compB) tcs_component_ids->elementData[j]--;
            }
        }

        if(other_pairs->size==0) free_alist_l(other_pairs);
        else{
            dist2pairs[M]=other_pairs;
            M--;
        }

        free_alist_l(vcptr);
    }
    free(dist2pairs);
    free_alist_i(tcs_component_ids);

    mylog("merging graph intermediate nodes...");
    int vertidx=seqCount;
    while(vertidx<all_verteres->size){
        Sequence *v=NULL;
        for(i=0,n=all_verteres->first;n;n=n->next,i++){
            if(i==vertidx){
                v=n->value;
                if((v->in_edges->size+v->out_edges->size)>2) vertidx++;
                else break;
            }
        }
        if((v->in_edges->size+v->out_edges->size)!=2){
            fprintf(stderr, "[error]Intermediate vertex has degree less than 2.\n");
            exit(0);
        }
        AList_l *edges=get_all_edges(v);
        Edge *e1=edges->elementData[0];
        Edge *e2=edges->elementData[1];
        free_alist_l(edges);
        Sequence *u=get_edge_opposite(e1, v);
        Sequence *w=get_edge_opposite(e2, v);
        if(u==w || v==w || u==v){
            fprintf(stderr, "[error]Unexpected multiple edges or self edge.\n");
            exit(0);
        }
        int weight=e1->weight+e2->weight;
        remove_vertex(v);
        if(u->index<w->index) add_edge(u, w, weight);
        else add_edge(w, u, weight);
    }
    for(i=0,n=all_verteres->first;n;n=n->next,i++){
        Sequence *v=n->value;
        v->index=i;
    }
}

AList_l *all_tcs_scores=NULL;

int findIntermediates(Pair *intPair, Sequence *u, Sequence *v, int dist){
    int maxScore=TCS_MIN_SCORE+1;
    int minPathLength=dist;

    int compU=tcs_component_ids->elementData[u->index];
    int compV=tcs_component_ids->elementData[v->index];

    if(compU==compV){fprintf(stderr, "[error]Attempting to find intermediates within a component.\n");exit(0);}

    intPair->left=u;
    intPair->right=v;

    int i, j;
    LListNode_l *n1, *n2;

    if(msn_thread_num>1) all_tcs_scores=new_alist_l(16);

    for(i=0,n1=all_verteres->first;n1&&i<tcs_component_ids->size;n1=n1->next,i++){
        if(tcs_component_ids->elementData[i]!=compU && tcs_component_ids->elementData[i]>=0)continue;
        int pathUI=get_path_length(u->index, i);
        if(pathUI==-1) continue;
        if(pathUI>=dist) continue;
        for(j=0,n2=all_verteres->first;n2&&j<tcs_component_ids->size;n2=n2->next,j++){
            if(tcs_component_ids->elementData[j]!=compV && tcs_component_ids->elementData[j]>=0)  continue;
            int pathVJ=get_path_length(v->index, j);
            if(pathVJ==-1) continue;
            if(pathVJ+pathUI>=dist) continue;
            int dP=dist-pathVJ-pathUI;
            if(msn_thread_num==1){
                int score=computeScore(((Sequence *)(n1->value))->index, ((Sequence *)(n2->value))->index, compU, compV, dP, dist);
                if(score>maxScore || (score==maxScore && dP<minPathLength)){
                    minPathLength = dP;
                    maxScore = score;
                    intPair->left=n1->value;
                    intPair->right=n2->value;
                }
            }else{
                TCS_tmp1 *tcs=my_new(1, sizeof(TCS_tmp1));
                tcs->index1=((Sequence *)(n1->value))->index;
                tcs->index2=((Sequence *)(n2->value))->index;
                tcs->compU=compU;
                tcs->compV=compV;
                tcs->dP=dP;
                tcs->clustDist=dist;
                tcs->left=n1->value;
                tcs->right=n2->value;
                alist_l_add(all_tcs_scores, tcs);
                thread_pool_add_worker(msn_thread_pool, computeScore_one, all_tcs_scores->size-1);
            }
        }
    }

    if(msn_thread_num>1){
        thread_pool_invoke_all(msn_thread_pool);
        for(i=0;i<all_tcs_scores->size;i++){
            TCS_tmp1 *tcs=all_tcs_scores->elementData[i];
            if(tcs->score>maxScore || (tcs->score==maxScore && tcs->dP<minPathLength)){
                minPathLength=tcs->dP;
                maxScore=tcs->score;
                intPair->left=tcs->left;
                intPair->right=tcs->right;
            }
            free(tcs);
        }
        free_alist_l(all_tcs_scores);
    }

    return minPathLength;
}

int multi_thread_pre;
int multi_thread_start_index;
int multi_thread_end_index;
int multi_thread_dist;

void newCompositePath(Sequence *start, Sequence *end, int dist){
    update_path_lengths_with_edge(start->index, end->index, dist);

    Sequence *u=start, *v;

    int i, j, pre=all_verteres->size;

    for(i=1;i<dist;i++){
        memset(loginfo, 0, (l_seq+1)*sizeof(char));
        sprintf(loginfo, "IN%d", msn_intermediate_index++);
        v=add_vertex(loginfo, loginfo, 1);
        v->index=all_verteres->size-1;
        alist_i_add(tcs_component_ids, -1);
        add_edge(u, v, 1);
        u=v;
    }
    add_edge(u, end, 1);

    ensure_path_lengths();
    if(msn_thread_num==1){
        for(j=0;j<pre;j++){
            int a=get_path_length(j, start->index);
            int b=get_path_length(j, end->index);
            for(i=pre;i<all_verteres->size;i++){
                int dis1=i-pre+1;
                int dis2=dist-dis1;
                if(a!=-1 && b!=-1) set_path_length(i, j, min(a+dis1, b+dis2));
            }
        }
    }else{
        multi_thread_pre=pre;
        multi_thread_start_index=start->index;
        multi_thread_end_index=end->index;
        multi_thread_dist=dist;
        split_span_for_threads(msn_thread_num, 0, pre);
        for(i=0;i<msn_thread_num;i++) thread_pool_add_worker(msn_thread_pool, newCompositePath_one, i);
        thread_pool_invoke_all(msn_thread_pool);
    }
    for(i=pre;i<all_verteres->size;i++){
        for(j=pre;j<=i;j++) set_path_length(i, j, i-j);
    }
}

void newCompositePath_one(int64_t para){
    int i, j, index=(int)(para&0xFFFFFFFFL);
    int start=starts_thread[index];
    int end=ends_thread[index];

    for(j=start;j<end;j++){
        int a=get_path_length(j, multi_thread_start_index);
        int b=get_path_length(j, multi_thread_end_index);
        for(i=multi_thread_pre;i<all_verteres->size;i++){
            int dis1=i-multi_thread_pre+1;
            int dis2=multi_thread_dist-dis1;
            if(a!=-1 && b!=-1) set_path_length(i, j, min(a+dis1, b+dis2));
        }
    }
}

void computeScore_one(int64_t para){
    TCS_tmp1 *tcs=all_tcs_scores->elementData[(int)(para&0xFFFFFFFFL)];
    tcs->score=computeScore(tcs->index1, tcs->index2, tcs->compU, tcs->compV, tcs->dP, tcs->clustDist);
}

int computeScore(int index1, int index2, int compU, int compV, int dP, int clustDist){
    int i, j, score=0;

    for(i=0;i<n_seq;i++){
        if(tcs_component_ids->elementData[i]!=compU) continue;
        for(j=0;j<n_seq;j++){
            if(tcs_component_ids->elementData[j]!=compV) continue;
            //--
            int a=get_path_length(index1, i);
            int b=get_path_length(index2, j);
            if(a==-1 || b==-1) score-=LONGPENALTY;
            else{
                int totalPath=dP+a+b;
                int dis=distances[(int64_t)i*(int64_t)r_distance+(int64_t)j];
                if (totalPath==dis) score+=BONUS;
                else if(totalPath>dis) score-=LONGPENALTY;
                else{
                    if(totalPath<clustDist) return TCS_MIN_SCORE;
                    else score-=SHORTCUTPENALTY;
                }
            }
        }
    }

    return score;
}

inline int get_path_length(int index1, int index2){
    return path_lengths[index1*r_path_lengths+index2];
}

inline void set_path_length(int index1, int index2, int dis){
    path_lengths[index1*r_path_lengths+index2]=path_lengths[index2*r_path_lengths+index1]=dis;
}

void ensure_path_lengths(){
    int seqCount=all_verteres->size;
    if(path_lengths==NULL){
        r_path_lengths=seqCount*5;
        path_lengths=my_new(r_path_lengths*r_path_lengths, sizeof(int));
        memset(path_lengths, -1, r_path_lengths*r_path_lengths*sizeof(int));
    }else if(r_path_lengths<seqCount){
        mylog("resizing...");
        int i, j, pre=r_path_lengths;
        r_path_lengths=seqCount*2;
        int *new_path_lengths=my_new(r_path_lengths*r_path_lengths, sizeof(int));
        memset(new_path_lengths, -1, r_path_lengths*r_path_lengths*sizeof(int));
        for(i=0;i<pre;i++){
            for(j=0;j<=i;j++) new_path_lengths[i*r_path_lengths+j]=new_path_lengths[j*r_path_lengths+i]=path_lengths[i*pre+j];
        }
        free(path_lengths);
        path_lengths=new_path_lengths;
    }
}

void update_path_lengths(){
    ensure_path_lengths();

    int i, j, k, seqCount=all_verteres->size;
    LListNode_l *n;

    for(i=0,n=all_verteres->first;n;n=n->next,i++){
        Sequence *v=n->value;
        v->index=i;
        set_path_length(i, i, 0);
    }

    for(n=all_edges->first;n;n=n->next){
        Edge *e=n->value;
        Sequence *u=e->from;
        Sequence *v=e->to;
        set_path_length(u->index, v->index, e->weight);
    }

    for(k=0;k<seqCount;k++){
        for(i=0;i<seqCount;i++){
            for(j=0;j<i;j++){
                int a=get_path_length(i, k);
                int b=get_path_length(k, j);
                if(a!=-1 && b!=-1){
                    int throughK=a+b;
                    int pre=get_path_length(i, j);
                    if(pre==-1 || throughK<pre) set_path_length(i, j, throughK);
                }
            }
        }
    }
}

int multi_thread_u_index;
int multi_thread_v_index;
int multi_thread_weight;
int multi_thread_k;
int multi_thread_finished_num;
AList_i *multi_thread_list0=NULL;
AList_i *multi_thread_list1=NULL;
AList_i *multi_thread_list2=NULL;

void update_path_lengths_with_edge(int u_index, int v_index, int weight){
    multi_thread_u_index=u_index;
    multi_thread_v_index=v_index;
    multi_thread_weight=weight;

    int dis=get_path_length(u_index, v_index);
    if(dis!=-1){
        set_path_length(u_index, v_index, min(dis, weight));

        if(multi_thread_list0==NULL) multi_thread_list0=new_alist_i(all_verteres->size);
        else multi_thread_list0->size=0;

        int i, j, k, l, m, n, t;

        for(i=0;i<all_verteres->size;i++){
            if(get_path_length(i, u_index)!=-1 || get_path_length(i, v_index)!=-1) alist_i_add(multi_thread_list0, i);
        }

        if(msn_thread_num==1){
            for(l=0;l<multi_thread_list0->size;l++){
                k=multi_thread_list0->elementData[l];
                for(m=0;m<multi_thread_list0->size;m++){
                    if(m==l) continue;
                    i=multi_thread_list0->elementData[m];
                    for(n=0;n<m;n++){
                        j=multi_thread_list0->elementData[n];
                        int a=get_path_length(i, k);
                        int b=get_path_length(k, j);
                        if(a!=-1 && b!=-1){
                            int throughK=a+b;
                            int pre=get_path_length(i, j);
                            if(pre==-1 || throughK<pre) set_path_length(i, j, throughK);
                        }
                    }
                }
            }
        }else{
            multi_thread_k=0;
            multi_thread_finished_num=0;
            for(t=0;t<msn_thread_num;t++) thread_pool_add_worker(msn_thread_pool, update_path_lengths_with_edge_one1, t);
            thread_pool_invoke_all(msn_thread_pool);
        }
        return;
    }

    int i, j;

    if(multi_thread_list1==NULL){
        multi_thread_list1=new_alist_i(16);
        multi_thread_list2=new_alist_i(16);
    }else{
        multi_thread_list1->size=0;
        multi_thread_list2->size=0;
    }

    for(i=0;i<all_verteres->size;i++){
        int dis1=get_path_length(i, u_index);
        int dis2=get_path_length(i, v_index);
        if(dis1!=-1 && dis2==-1) alist_i_add(multi_thread_list1, i);
        if(dis1==-1 && dis2!=-1) alist_i_add(multi_thread_list2, i);
    }

    if(msn_thread_num==1){
        for(i=0;i<multi_thread_list1->size;i++){
            int index1=multi_thread_list1->elementData[i];
            int dis1=get_path_length(index1, u_index)+weight;
            for(j=0;j<multi_thread_list2->size;j++){
                int index2=multi_thread_list2->elementData[j];
                set_path_length(index1, index2, dis1+get_path_length(index2, v_index));
            }
        }
    }else{
        if(multi_thread_list1->size<multi_thread_list2->size){
            AList_i *tmp1=multi_thread_list1;
            multi_thread_list1=multi_thread_list2;
            multi_thread_list2=tmp1;
            int tmp2=multi_thread_u_index;
            multi_thread_u_index=multi_thread_v_index;
            multi_thread_v_index=tmp2;
        }
        split_span_for_threads(msn_thread_num, 0, multi_thread_list1->size);
        for(i=0;i<msn_thread_num;i++) thread_pool_add_worker(msn_thread_pool, update_path_lengths_with_edge_one2, i);
        thread_pool_invoke_all(msn_thread_pool);
    }
}

void update_path_lengths_with_edge_one1(int64_t para){
    int i, j, k, m, n, index=(int)(para&0xFFFFFFFFL);
    while(1){
        k=multi_thread_list0->elementData[multi_thread_k];
        for(m=multi_thread_list0->size-index-1;m>=0;m-=msn_thread_num){
            if(m==multi_thread_k) continue;
            i=multi_thread_list0->elementData[m];
            for(n=0;n<m;n++){
                j=multi_thread_list0->elementData[n];
                int a=get_path_length(i, k);
                int b=get_path_length(k, j);
                if(a!=-1 && b!=-1){
                    int throughK=a+b;
                    int pre=get_path_length(i, j);
                    if(pre==-1 || throughK<pre) set_path_length(i, j, throughK);
                }
            }
        }
        pthread_mutex_lock(lock_for_multi_thread);
        multi_thread_finished_num++;
        if(multi_thread_finished_num<msn_thread_num) pthread_cond_wait(cond_for_multi_thread, lock_for_multi_thread);
        else{
            multi_thread_k++;
            multi_thread_finished_num=0;
            pthread_cond_broadcast(cond_for_multi_thread);
        }
        pthread_mutex_unlock(lock_for_multi_thread);
        if(multi_thread_k==multi_thread_list0->size) break;
    }
}

void update_path_lengths_with_edge_one2(int64_t para){
    int i, j, index=(int)(para&0xFFFFFFFFL);
    int start=starts_thread[index];
    int end=ends_thread[index];

    for(i=start;i<end;i++){
        int index1=multi_thread_list1->elementData[i];
        int dis1=get_path_length(index1, multi_thread_u_index)+multi_thread_weight;
        for(j=0;j<multi_thread_list2->size;j++){
            int index2=multi_thread_list2->elementData[j];
            set_path_length(index1, index2, dis1+get_path_length(index2, multi_thread_v_index));
        }
    }
}

//=============================================

void compute_msn_only(){
    sprintf(loginfo, "----------computing %d MSN...", msn_compute_number++);mylog(loginfo);

    int _strict=epsilon==0 ? 1:0;

    int i, j, k;
    LListNode_l *n, *n1, *n2;

    calculate_distances();

    int seqCount=all_verteres->size;
    int ncomps=all_verteres->size;

    int *component=my_new(seqCount, sizeof(int));
    int64_t *dist2pairs=my_new(max_distance+1, sizeof(int64_t));

    int max_value=1<<30;

    for(i=0,n1=all_verteres->first;n1;i++,n1=n1->next){
        component[i]=i;
        Sequence *v1=n1->value;
        for(j=0,n2=all_verteres->first;n2&&j<i;j++,n2=n2->next){
            Sequence *v2=n2->value;
            int dis=distances[(int64_t)i*(int64_t)r_distance+(int64_t)j];
            AList_l *list=dist2pairs[dis];
            if(list==NULL){
                list=new_alist_l(16);
                dist2pairs[dis]=list;
            }
            Pair *p=my_new(1, sizeof(Pair));
            p->left=v1;
            p->right=v2;
            alist_l_add(list, p);
        }
    }

    int threshold;
    for(threshold=0;threshold<=max_distance;threshold++){
        AList_l *vcptr=dist2pairs[threshold];
        if(vcptr==NULL) continue;
        if(threshold>max_value) break;
        //--
        AList_l *new_pairs=new_alist_l(vcptr->size);
        for(i=0;i<vcptr->size;i++){
            Pair *pair=vcptr->elementData[i];
            Sequence *u=pair->left;
            Sequence *v=pair->right;
            if(!_strict || (component[u->index]!=component[v->index])){
                Edge *e=add_edge(u, v, threshold);
                alist_l_add(new_pairs, pair);
            }else{
                free(pair);
            }
        }
        //--
        for(i=0;i<new_pairs->size;i++){
            Pair *pair=new_pairs->elementData[i];
            Sequence *u=pair->left;
            Sequence *v=pair->right;
            int compU=component[u->index];
            int compV=component[v->index];
            if(compU!=compV){
                if(compU<compV){
                    int tmp=compU;
                    compU=compV;
                    compV=tmp;
                }
                for(j=0;j<seqCount;j++){
                    if(component[j]==compU) component[j]=compV;
                    else if(component[j]>compU) component[j]--;
                }
                ncomps--;
            }
            if(ncomps==1 && max_value==(1<<30)) max_value=threshold+epsilon;
            free(pair);
        }
        free_alist_l(new_pairs);
        free_alist_l(vcptr);
        dist2pairs[threshold]=NULL;
    }

    free(dist2pairs);
    free(component);
}

//=============================================

int mjn_min_cost;
int *mjn_min_cost_array;

void compute_mjn_graph(){
    compute_mjn();

    int i, changed=0;
    LListNode_l *n;
    do{
        AList_l *feasibleLinks=new_alist_l(16);
        compute_msn(feasibleLinks);

        clear_all_edges();

        for(i=0;i<feasibleLinks->size;i++){
            Edge *e=feasibleLinks->elementData[i];
            add_edge(e->from, e->to, e->weight);
            free(e);
        }
        free_alist_l(feasibleLinks);

        changed=remove_obsolete_vertexes();
  }while(changed);
}

void compute_msn(AList_l *feasibleLinks){
    sprintf(loginfo, "----------calculating %d MSN...", msn_compute_number++);mylog(loginfo);

    int i, j, k, t;
    LListNode_l *n, *n1, *n2;

    calculate_distances();

    //-- clear all edges
    clear_all_edges();

    int seqCount=all_verteres->size;
    int ncomps=all_verteres->size;

    int avx_count=seqCount/16;
    if((avx_count*16)<seqCount) avx_count++;

    uint16_t *msnComp=calloc_align_ptr(avx_count*16*sizeof(uint16_t), 32);
    uint16_t *thresholdComp=calloc_align_ptr(avx_count*16*sizeof(uint16_t), 32);
    int64_t *dist2pairs=my_new(max_distance+1, sizeof(int64_t));

    __m256i _all=_mm256_set1_epi16(0xFFFF);

    int max_value=1<<30;

    for(i=0,n1=all_verteres->first;n1;i++,n1=n1->next){
        msnComp[i]=i;
        thresholdComp[i]=i;
        Sequence *v1=n1->value;
        for(j=0,n2=all_verteres->first;n2&&j<i;j++,n2=n2->next){
            Sequence *v2=n2->value;
            int dis=distances[(int64_t)i*(int64_t)r_distance+(int64_t)j];
            AList_l *list=dist2pairs[dis];
            if(list==NULL){
                list=new_alist_l(16);
                dist2pairs[dis]=list;
            }
            Pair *p=my_new(1, sizeof(Pair));
            p->left=v1;
            p->right=v2;
            alist_l_add(list, p);
        }
    }

    int threshold;
    for(threshold=min_distance;threshold<=max_distance;threshold++){
        AList_l *vcptr=dist2pairs[threshold];
        if(vcptr==NULL) continue;
        if(threshold>max_value) break;
        //--
        for(i=0;i<seqCount;i++){
            for(j=0;j<i;j++){
                if(thresholdComp[i]!=thresholdComp[j] && distances[(int64_t)i*(int64_t)r_distance+(int64_t)j]<(threshold-epsilon)){
                    uint16_t oldComp = thresholdComp[i];
                    uint16_t newComp = thresholdComp[j];
                    if(oldComp<newComp){
                        oldComp=thresholdComp[j];
                        newComp=thresholdComp[i];
                    }
                    if(is_use_avx){
                        __m256i _old=_mm256_set1_epi16(oldComp);
                        __m256i _new=_mm256_set1_epi16(newComp);
                        uint16_t *copy=thresholdComp;
                        for(k=0;k<avx_count;k++){
                            __m256i a=_mm256_load_si256(copy);
                            __m256i b=_mm256_cmpeq_epi16(a, _old);
                            _mm256_store_si256(copy, _mm256_or_si256(_mm256_and_si256(b, _new), _mm256_and_si256(_mm256_xor_si256(b, _all), a)));
                            copy+=16;
                        }
                    }else{
                        for(k=0;k<seqCount;k++){
                            if(thresholdComp[k]==oldComp) thresholdComp[k]=newComp;
                        }
                    }
                }
            }
        }
        //--
        for(i=0;i<vcptr->size;i++){
            Pair *pair=vcptr->elementData[i];
            Sequence *u=pair->left;
            Sequence *v=pair->right;
            if(mjn_compute_number>1){   // is called by mjn
                if(thresholdComp[u->index]!=thresholdComp[v->index]){
                    Edge *e=my_new(1, sizeof(Edge));
                    e->from=u;
                    e->to=v;
                    e->weight=threshold;
                    alist_l_add(feasibleLinks, e);
                }
            }else{
                Edge *e=add_edge(u, v, threshold);
                if(thresholdComp[u->index]!=thresholdComp[v->index]) alist_l_add(feasibleLinks, e);
            }
        }
        //--
        for(i=0;i<vcptr->size;i++){
            Pair *pair=vcptr->elementData[i];
            Sequence *u=pair->left;
            Sequence *v=pair->right;
            uint16_t compU=msnComp[u->index];
            uint16_t compV=msnComp[v->index];
            if(compU!=compV){
                if(compU<compV){
                    uint16_t tmp=compU;
                    compU=compV;
                    compV=tmp;
                }
                if(is_use_avx){
                    __m256i _old=_mm256_set1_epi16(compU);
                    __m256i _new=_mm256_set1_epi16(compV);
                    uint16_t *copy=msnComp;
                    for(j=0;j<avx_count;j++){
                        __m256i a=_mm256_load_si256(copy);
                        __m256i b=_mm256_cmpeq_epi16(a, _old);
                        _mm256_store_si256(copy, _mm256_or_si256(_mm256_and_si256(b, _new), _mm256_and_si256(_mm256_xor_si256(b, _all), a)));
                        copy+=16;
                    }
                }else{
                    for(j=0;j<seqCount;j++){
                        if(msnComp[j]==compU) msnComp[j]=compV;
                    }
                }
                ncomps--;
            }
            if(ncomps==1 && max_value==(1<<30)) max_value=threshold+epsilon;
            free(pair);
        }
        //--
        free_alist_l(vcptr);
        dist2pairs[threshold]=NULL;
    }

    free(dist2pairs);
    free_align_ptr(msnComp);
    free_align_ptr(thresholdComp);

    sprintf(loginfo, "threshold=%d\tfeasibleLinks=%d", threshold, feasibleLinks->size);mylog(loginfo);
}

void compute_mjn(){
    sprintf(loginfo, "==========calculating %d MJN...", mjn_compute_number++);mylog(loginfo);
    msn_compute_number=1;

    int changed;

    int i, j, k, l, x, y;
    LListNode_l *n, *n1, *n2, *n3;

    do{
        AList_l *feasibleLinks=new_alist_l(16);
        compute_msn(feasibleLinks);

        clear_all_edges();

        AList_l *tmp_feasibleLinks=new_alist_l(feasibleLinks->size);
        for(i=0;i<feasibleLinks->size;i++){
            Edge *e=feasibleLinks->elementData[i];
            alist_l_add(tmp_feasibleLinks, add_edge(e->from, e->to, e->weight));
            free(e);
        }
        free_alist_l(feasibleLinks);
        feasibleLinks=tmp_feasibleLinks;

        changed=remove_obsolete_vertexes();

        mylog("calculating med vertex...");

        all_meds=new_alist_l(16);
        l_table_medes=feasibleLinks->size*100;
        table_medes=my_new(l_table_medes, sizeof(int64_t));

        for(i=0;i<feasibleLinks->size;i++){
            Edge *e1=feasibleLinks->elementData[i];
            Sequence *u=e1->from;
            Sequence *v=e1->to;
            for(j=0;j<i;j++){
                Edge *e2=feasibleLinks->elementData[j];
                Sequence *w;
                if(e2->from==u || e2->from==v) w=e2->to;
                else if(e2->to==u || e2->to==v) w=e2->from;
                else continue;
                Entry_med *entry=put_hash_meds(table_medes, l_table_medes, u, v, w, NULL, NULL);
                if(entry->freq==1) alist_l_add(all_meds, entry);
            }
        }
        free_alist_l(feasibleLinks);
        sprintf(loginfo, "all_med_triple=%d", all_meds->size);mylog(loginfo);

        mjn_min_cost=1<<30;
        if(1){
            split_span_for_threads(msn_thread_num, 0, all_meds->size);
            for(i=0;i<msn_thread_num;i++) thread_pool_add_worker(msn_thread_pool, calculate_min_cost_one, i);
            thread_pool_invoke_all(msn_thread_pool);
        }else{
            mjn_min_cost_array=my_new(msn_thread_num, sizeof(int));
            for(i=0;i<msn_thread_num;i++) mjn_min_cost_array[i]=1<<30;
            for(i=0;i<all_meds->size;i++) thread_pool_add_worker(msn_thread_pool, calculate_min_cost_one2, i);
            thread_pool_invoke_all(msn_thread_pool);
            for(i=0;i<msn_thread_num;i++){
                if(mjn_min_cost>mjn_min_cost_array[i]) mjn_min_cost=mjn_min_cost_array[i];
            }
            free(mjn_min_cost_array);
        }

        sprintf(loginfo, "min_cost=%d", mjn_min_cost);mylog(loginfo);

        int add_vertex_number=0;
        for(i=0;i<all_meds->size;i++){
            Entry_med *entry=all_meds->elementData[i];
            AList_l *medes=entry->medes;
            AList_i *costs=entry->costs;
            for(j=0;j<medes->size;j++){
                char *str=medes->elementData[j];
                int cost=costs->elementData[j];
                if(cost<=(mjn_min_cost+epsilon)){
                    if(get_hash_seqs(table_verteres, l_table_verteres, str, l_seq)==NULL){
                        sprintf(loginfo, "IN%d", msn_intermediate_index++);
                        add_vertex(loginfo, str, 1);
                        add_vertex_number++;
                        changed=1;
                    }
                }
                free(str);
            }
            free_alist_l(medes);
            free_alist_i(costs);
            free(entry);
        }

        free_alist_l(all_meds);
        free(table_medes);

        sprintf(loginfo, "add med vertex: %d", add_vertex_number);mylog(loginfo);
    }while(changed);
}

void calculate_min_cost_one(int64_t para){
    int thread_index=(int)(para&0xFFFFFFFFL);
    int start=starts_thread[thread_index];
    int end=ends_thread[thread_index];

    int min_cost=1<<30;

    int i, j;
    for(i=start;i<end;i++){
        Entry_med *entry=all_meds->elementData[i];
        Sequence *u=entry->key1;
        Sequence *v=entry->key2;
        Sequence *w=entry->key3;
        if(u==v || u==w || v==w){fprintf(stderr, "[error]here1!\n");exit(0);}
        //--
        AList_l *medes=new_alist_l(16);
        AList_i *costs=new_alist_i(16);
        SBuilder *sb=new_s_builder(l_seq);
        computeQuasiMedianSeqs(u->seq, v->seq, w->seq, sb, l_seq, medes, costs);
        free_s_builder(sb);
        //--
        for(j=0;j<costs->size;j++){
            char *str=medes->elementData[j];
            if(get_hash_seqs(table_verteres, l_table_verteres, str, l_seq)==NULL){
                int cost=costs->elementData[j];
                if(min_cost>cost) min_cost=cost;
            }
        }
        //--
        entry->medes=medes;
        entry->costs=costs;
    }

    pthread_mutex_lock(lock_for_multi_thread);
    if(mjn_min_cost>min_cost) mjn_min_cost=min_cost;
    pthread_mutex_unlock(lock_for_multi_thread);
}

void calculate_min_cost_one2(int64_t para){
    int thread_index=(int)(para>>32);
    int i=(int)(para&0xFFFFFFFFL);

    int min_cost=1<<30;

    int j;
    {
        Entry_med *entry=all_meds->elementData[i];
        Sequence *u=entry->key1;
        Sequence *v=entry->key2;
        Sequence *w=entry->key3;
        if(u==v || u==w || v==w){fprintf(stderr, "[error]here1!\n");exit(0);}
        //--
        AList_l *medes=new_alist_l(16);
        AList_i *costs=new_alist_i(16);
        SBuilder *sb=new_s_builder(l_seq);
        computeQuasiMedianSeqs(u->seq, v->seq, w->seq, sb, l_seq, medes, costs);
        free_s_builder(sb);
        //--
        for(j=0;j<costs->size;j++){
            char *str=medes->elementData[j];
            if(get_hash_seqs(table_verteres, l_table_verteres, str, l_seq)==NULL){
                int cost=costs->elementData[j];
                if(min_cost>cost) min_cost=cost;
            }
        }
        //--
        entry->medes=medes;
        entry->costs=costs;
    }

    if(mjn_min_cost_array[thread_index]>min_cost) mjn_min_cost_array[thread_index]=min_cost;
}

void computeQuasiMedianSeqs(char *seqA, char *seqB, char *seqC, SBuilder *sb, int len, AList_l *list, AList_i *costs){
    if(sb->size==len){
        alist_l_add(list, str_copy_with_len(sb->str, sb->size));
        alist_i_add(costs, computeCost(seqA, seqB, seqC, sb->str));
        return;
    }
    char c1=seqA[sb->size];
    char c2=seqB[sb->size];
    char c3=seqC[sb->size];
    if(c1==c2 || c1==c3){
        s_builder_add_char(sb, c1);
        computeQuasiMedianSeqs(seqA, seqB, seqC, sb, len, list, costs);
        sb->size--;
    }else if(c2==c3){
        s_builder_add_char(sb, c2);
        computeQuasiMedianSeqs(seqA, seqB, seqC, sb, len, list, costs);
        sb->size--;
    }else{
        s_builder_add_char(sb, c1);
        computeQuasiMedianSeqs(seqA, seqB, seqC, sb, len, list, costs);
        sb->size--;
        s_builder_add_char(sb, c2);
        computeQuasiMedianSeqs(seqA, seqB, seqC, sb, len, list, costs);
        sb->size--;
        s_builder_add_char(sb, c3);
        computeQuasiMedianSeqs(seqA, seqB, seqC, sb, len, list, costs);
        sb->size--;
    }
}

inline int computeCost(char *seqA, char *seqB, char *seqC, char *med){
    return pairwiseDistance(seqA, med)+pairwiseDistance(seqB, med)+pairwiseDistance(seqC, med);
}

int pairwiseDistance(char *s1, char *s2){
    int k, num=0;
    for(k=0;k<l_seq;k++){
        char c1=s1[k];
        char c2=s2[k];
        if(c1=='?' || c2=='?') continue;
        if(c1!=c2) num+=weigths[k];
    }
    return num;
}

int remove_obsolete_vertexes(){
    sprintf(loginfo, "removing obsolete vertexes...");mylog(loginfo);

    int changed=1;
    int vertsRemoved=0;
    LListNode_l *n;

    int i, j;

    while(changed){
        AList_l *obsoleteVerts=new_alist_l(16);
        for(i=0,n=all_verteres->first;n;n=n->next,i++){
            if(i<n_seq) continue;
            Sequence *v=n->value;
            if((v->in_edges->size+v->out_edges->size)<2) alist_l_add(obsoleteVerts, v);
        }
        //--
        sprintf(loginfo, "removed_obsolete_vertexes=%d", obsoleteVerts->size);mylog(loginfo);
        //--
        if(obsoleteVerts->size==0) changed=0;
        else{
            vertsRemoved=1;
            for(i=0;i<obsoleteVerts->size;i++){
                Sequence *v=obsoleteVerts->elementData[i];
                remove_vertex(v);
            }
        }
        free_alist_l(obsoleteVerts);
    }

    return vertsRemoved;
}

//=============================================

void msn_output_graph(){
    int k, len;

    sprintf(loginfo, "%s.gml", msn_outFile);
    GzStream *out=gz_stream_open(loginfo, "w");

    gz_write(out, "graph [\n", 8);
    gz_write(out, "   directed 0\n", 14);

    LListNode_l *n;
    for(n=all_verteres->first;n;n=n->next){
        Sequence *v=n->value;
        //len=sprintf(loginfo, "vertex:\t%d\t%s\t%d\t%d\n", v->index, v->name, v->index<n_seq?0:1, v->dup_num);
        //gz_write(out2, loginfo, len);
        gz_write(out, "   node [\n", 10);
        gz_write(out, "      data [\n", 13);
        gz_write(out, loginfo, sprintf(loginfo, "         Frequency \"frequency=%.1f\n", v->index>=n_seq ? 1.0:(double)(v->dup_num)));
        for(k=0;k<v->dupnames->size;k++) gz_write(out, loginfo, sprintf(loginfo, "%s\n", v->dupnames->elementData[k]));
        gz_write(out, "\"\n      ]\n", 10);
        gz_write(out, loginfo, sprintf(loginfo, "      id %d\n   ]\n", v->index));
    }

    for(n=all_edges->first;n;n=n->next){
        Edge *edge=n->value;
        //len=sprintf(loginfo, "edge:\t%d\t%d\t%d\n", edge->from->index, edge->to->index, edge->weight);
        //gz_write(out2, loginfo, len);
        gz_write(out, "   edge [\n", 10);
        gz_write(out, "      data [\n", 13);
        gz_write(out, loginfo, sprintf(loginfo, "         Changes \"%d\"\n", edge->weight));
        gz_write(out, "      ]\n", 8);
        gz_write(out, loginfo, sprintf(loginfo, "      source %d\n", edge->from->index));
        gz_write(out, loginfo, sprintf(loginfo, "      target %d\n", edge->to->index));
        gz_write(out, "   ]\n", 5);
    }

    gz_write(out, "]\n", 2);

    gz_stream_destory(out);

    sprintf(loginfo, "%s.json", msn_outFile);
    out=gz_stream_open(loginfo, "w");

    gz_write(out, "{\n", 2);

    gz_write(out, "  \"nodes\": [\n", 13);
    int i;
    for(i=0,n=all_verteres->first;n;n=n->next,i++){
        Sequence *v=n->value;
        SBuilder *sb=new_s_builder(16);
        for(k=0;k<v->dupnames->size;k++){
            s_builder_add_str(sb, v->dupnames->elementData[k]);
            if(k<(v->dupnames->size-1)) s_builder_add_char(sb, ';');
        }
        double freq=v->index>=n_seq ? 1.0:(double)(v->dup_num);
        gz_write(out, loginfo, sprintf(loginfo, "    {\"id\": %d, \"frequency\": %f, \"title1\": \"%s\", \"title2\": \"%s\"}", v->index, freq, v->dupnames->elementData[0], sb->str));
        if(i<(all_verteres->size-1)) gz_write_char(out, ',');
        gz_write_char(out, '\n');
        free_s_builder(sb);
    }
    gz_write(out, "  ],\n", 5);

    gz_write(out, "  \"edges\": [\n", 13);
    for(i=0,n=all_edges->first;n;n=n->next,i++){
        Edge *edge=n->value;
        gz_write(out, loginfo, sprintf(loginfo, "    {\"source\": %d, \"target\": %d, \"change\": %d}", edge->from->index, edge->to->index, edge->weight));
        if(i<(all_edges->size-1)) gz_write_char(out, ',');
        gz_write_char(out, '\n');
    }
    gz_write(out, "  ]\n", 4);

    gz_write(out, "}\n", 2);

    gz_stream_destory(out);
}

//=============================================


