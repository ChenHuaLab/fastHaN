#include "tcs.h"

//=================================

ThreadPool *threadPool;
AList_l *allvect_for_multi_thread;
int tabs[1000];
int graphExists=0;
int numHapInt=0;
int numHapHap=0;
int internalNodeNumber=1;

//=================================

int main_tcs(int argc, char **argv){
    load_parameters(argc, argv);

    mylog("starting...");

    threadPool=new_thread_pool(thread_num);
    locker=my_new(1, sizeof(pthread_mutex_t));
    pthread_mutex_init(locker, NULL);
    evaluateMetric_results=my_new(thread_num, sizeof(int));
    edges=new_alist_l(16);

    tcs_read_input_file();
    calculate_connection_limit();
    if(components->size<2){
        fprintf(stdout, "the seq(marked duplication)'s number is only one!");
        return 0;
    }
    build_matrix();

    int i, j, graphNumber=0;

    mylog("connecting taxa...");
    for(i=0;i<realtaxa->size;i++){
        connectTaxa(realtaxa->elementData[i]);
        graphNumber++;
    }

    mylog("connecting components...");
    while(components->size>1){
        graphNumber++;
        if(!connectComponents2()) break;
    }

    if(graphExists){
        if(!is_merge_intermediate){
            mylog("calculating out group weights...");
            calculateOutgroupWeights();
            mylog("outputting graphes...");
            printGraph();
        }else{
            output_graph2();
        }
    }

    mylog("finished!");
}

void load_parameters(int argc, char **argv){
    argc--;
    argv++;

    if(argc<4 || argc>10) print_usage();
    if(argc%2==1) print_usage();

    thread_num=8;
    is_mask_ambiguous=0;
    is_merge_intermediate=0;

    int i;
    for(i=0;i<argc;i++){
        if(strcmp(argv[i], "-i")==0){
            i++;
            inFile=argv[i];
        }else if(strcmp(argv[i], "-t")==0){
            i++;
            thread_num=atoi(argv[i]);
        }else if(strcmp(argv[i], "-a")==0){
            i++;
            is_mask_ambiguous=atoi(argv[i]);
        }else if(strcmp(argv[i], "-m")==0){
            i++;
            is_merge_intermediate=atoi(argv[i]);
        }else if(strcmp(argv[i], "-o")==0){
            i++;
            outFile=argv[i];
        }else{
            fprintf(stderr, "'%s' parameter is not exist, please check!\n");
            exit(0);
        }
    }

    if(thread_num<1) thread_num=1;
}

void print_usage(){
    printf("===========================\n");
    printf("Author:   ChiLianjiang\n");
    printf("E-mail:   chilianjiang@126.com\n");
    printf("Date:     2021-03-24\n");
    printf("Version:  1.0\n\n");
    printf("Usage: fastHapNetwork original_tcs [arguments]\n");
    printf("  input:\n");
    printf("    -i    input phylip format file\n");
    printf("  options:\n");
    printf("    -t    (int)thread number(default:8)\n");
    printf("    -a    (int)is mark the sites which contains ambiguous base, 1:mask, 0:not(default:0)\n");
    printf("    -m    (int)is merge intermediate vertex, 1:merge, 0:not:(default:0)\n");
    printf("  output:\n");
    printf("    -o    output prefix(graph file)(GML and json format)\n\n");
    exit(0);
}

//=================================

void tcs_read_input_file(){
    sprintf(loginfo, "reading (%s)...", inFile);mylog(loginfo);

    int i, j, t, cnt=0, IUPACcounter=0;

    components=new_alist_l(16);
    alltaxa=new_alist_l(16);
    realtaxa=new_alist_l(16);

    AList_l *names=new_alist_l(16);
    AList_l *seqs=new_alist_l(16);
    tcs_get_data(inFile, names, seqs);

    for(t=0;t<names->size;t++){
        char *name=names->elementData[t];
        char *seq=seqs->elementData[t];
        //--
        TaxaItem *temptaxa=new_taxaitem2(name, seq_len, cnt);
        temptaxa->resolved=1;
        for(i=0;i<seq_len;i++){
            char c=toupper(seq[i]);
            if (c == 'R' || c == 'M' || c == 'W' || c == 'S' || c == 'K' || c == 'Y' ||
						c == 'H' || c == 'V' || c == 'D' || c == 'B' || c == 'X' || c == 'N') {
                c='?';
                if (IUPACcounter == 0) IUPACcounter++;
            }
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != '-' && c != '.' && c != '?'){
                fprintf(stderr, "error base=%d(c=%c), name=%s\n", c, c, name);
                exit(0);
            }
            temptaxa->characters[i]=c;
        }
        //--
        int is_dup=0;
        char newc, oldc;
        Component *tempc=NULL;
        TaxaItem *firsttaxa=NULL;
        if(components->size>0){
            Component *com=components->elementData[0];
            firsttaxa=com->taxa->elementData[0];
        }
        for(i=0;i<components->size;i++){
            tempc=components->elementData[i];
            int dup=1;
            TaxaItem *tempt=tempc->taxa->elementData[0];
            for(j=0;j<seq_len;j++) {
                newc = temptaxa->characters[j];
                oldc = tempt->characters[j];
                if (newc == '.') newc = firsttaxa->characters[j];
                if (oldc == '.') oldc = firsttaxa->characters[j];
                if (newc == '?' || oldc == '?') continue;
                if (newc != oldc) {dup = 0;break;}
            }
            if(dup){
                tempt->numduplicates++;
                is_dup=1;
                alist_l_add(tempt->dupnames, str_copy(name));
                free_taxaitem(temptaxa);
                break;
            }
        }
        if(!is_dup){
            tempc=new_component(cnt);
            alist_l_add(tempc->taxa, temptaxa);
            temptaxa->parentComponent=tempc;
            alist_l_add(alltaxa, temptaxa);
            alist_l_add(realtaxa, temptaxa);
            cnt++;
            alist_l_add(components, tempc);
        }
        //--
        free(name);
        free(seq);
    }
    free_alist_l(names);
    free_alist_l(seqs);

    sprintf(loginfo, "seq_num=%d\tcomponents_size=%d", seq_num, components->size);mylog(loginfo);
}

void tcs_get_data(char *file, AList_l *names, AList_l *seqs){
    int i, j, k, t, line_len, max_len=1000000000;

    GzStream *gz1=gz_stream_open(inFile, "r");

    //-- first line
    gz_read_util(gz1, '\n', line, max_len, &line_len);
    line_len=chmop_with_len(line, line_len);
    str_replace_char_with_no_copy(line, ' ', '\t');
    i=0;
    while(line[++i]!='\t');
    line[i]='\0';
    seq_num=atoi(line);
    seq_len=atoi(line+i+1);

    int new_seq_len=0;
    int *mask=tcs_get_data_mask(file, seq_len);

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
        alist_l_add(names, str_copy_with_len(name, l_name));
        if(!is_mask_ambiguous) alist_l_add(seqs, str_copy_with_len(seq, seq_len));
        else{
            SBuilder *sb=new_s_builder(seq_len);
            for(i=0;i<seq_len;i++){
                if(!mask[i]) s_builder_add_char(sb, toupper(seq[i]));
            }
            alist_l_add(seqs, sb->str);
            new_seq_len=sb->size;
            free(sb);
        }
    }
    gz_stream_destory(gz1);
    free(mask);

    if(!is_mask_ambiguous) return;

    sprintf(loginfo, "seq_len=%d\tmark_len=%d", seq_len, new_seq_len);mylog(loginfo);
    seq_len=new_seq_len;
}

int *tcs_get_data_mask(char *file, int l_mask){
    int *mask=my_new(l_mask, sizeof(int));

    int i, line_len, max_len=1000000000;

    tcs_is_nucleotide_base=1;

    GzStream *gz1=gz_stream_open(file, "r");
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
        for(i=0;i<l_mask;i++){
            char c=toupper(seq[i]);
            if (c == '-' || c == '.' ||
                c == 'R' || c == 'M' || c == 'W' || c == 'S' || c == 'K' || c == 'Y' ||
                c == 'H' || c == 'V' || c == 'D' || c == 'B' || c == 'X' || c == 'N'){
                c='?';
            }
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != '?'){
                tcs_is_nucleotide_base=0;
                goto pos;
            }
        }
    }
    pos:
    gz_stream_destory(gz1);

    gz1=gz_stream_open(file, "r");
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
        for(i=0;i<l_mask;i++){
            char c=toupper(seq[i]);
            if(tcs_is_nucleotide_base){
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

//=================================

MyIterator *new_my_iterator(AList_l *list){
    MyIterator *iter=my_new(1, sizeof(MyIterator));
    iter->count=0;
    iter->list=list;
    return iter;
}

int my_iterator_has_more(MyIterator *iter){
    return iter->count<iter->list->size ? 1:0;
}

void *my_iterator_next(MyIterator *iter){
    if(iter->count<iter->list->size){
        void *res=iter->list->elementData[iter->count++];
        return res;
    }
    return NULL;
}

TaxaItem *new_taxaitem1(char *n, int ident, char *chars, int l_chars){
    TaxaItem *item=my_new(1, sizeof(TaxaItem));

    item->paths=new_alist_l(16);
    item->name=str_copy(n);
    if(chars) item->characters=str_copy_with_len(chars, l_chars);
    item->l_characters=l_chars;
    item->id=ident;
    item->levelNumber=-1;

    item->realdist=new_alist_l(16);
    item->compdist=new_alist_l(16);
    item->nbor=new_alist_l(16);
    item->newconnections=new_alist_l(16);
    item->dupnames=new_alist_l(16);
    item->metricdist=new_alist_l(16);

    return item;
}

TaxaItem *new_taxaitem2(char *n, int length, int ident){
    TaxaItem *item=my_new(1, sizeof(TaxaItem));

    item->paths=new_alist_l(16);
    item->name=str_copy(n);
    item->characters=my_new(length, sizeof(char));
    item->l_characters=length;
    item->id=ident;
    item->levelNumber=-1;

    item->realdist=new_alist_l(16);
    item->compdist=new_alist_l(16);
    item->nbor=new_alist_l(16);
    item->newconnections=new_alist_l(16);
    item->dupnames=new_alist_l(16);
    item->metricdist=new_alist_l(16);

    return item;
}

void free_taxaitem(TaxaItem *item){
    int i, j;
    if(item->paths){

        free_alist_l(item->paths);
    }
    if(item->name) free(item->name);
    if(item->characters) free(item->characters);
    if(item->realdist){

        free_alist_l(item->realdist);
    }
    if(item->compdist){

        free_alist_l(item->compdist);
    }
    if(item->nbor){

        free_alist_l(item->nbor);
    }
    if(item->newconnections){

        free_alist_l(item->newconnections);
    }
    if(item->dupnames){
        for(i=0;i<item->dupnames->size;i++){
            char *name=item->dupnames->elementData[i];
            free(name);
        }
        free_alist_l(item->dupnames);
    }
    if(item->metricdist){

        free_alist_l(item->metricdist);
    }
    free(item);
}

Component *new_component(int ident){
    Component *com=my_new(1, sizeof(Component));
    com->id=ident;
    com->taxa=new_alist_l(16);
    com->compdist=new_alist_l(16);
    com->mindist=new_distance1();
    return com;
}

Distance *new_distance1(){
    return my_new(1, sizeof(Distance));
}

Distance *new_distance2(int src, int dest, int dist){
    Distance *dis=my_new(1, sizeof(Distance));
    dis->source=src;
    dis->destination=dest;
    dis->distance=dist;
    return dis;
}

Distance *new_distance3(int src, int dest, int dist, Component *srcc, Component *dstc){
    Distance *dis=my_new(1, sizeof(Distance));
    dis->source=src;
    dis->destination=dest;
    dis->distance=dist;
    dis->sc=srcc;
    dis->dc=dstc;
    return dis;
}

inline void clone_distance(Distance *dest, Distance *src){
    memcpy(dest, src, sizeof(Distance));
}

Path *new_path(TaxaItem *mySourcet, TaxaItem *myDestt, int myType){
    Path *path=my_new(1, sizeof(Path));
    path->source=mySourcet;
    path->dest=myDestt;
    path->type=myType;
    path->resolved=0;
    path->edges=new_alist_l(16);
    path->differences=new_alist_l(16);
    path->isAmbiguous=0;
    return path;
}

LWEdge *new_LWEdge(TaxaItem *myDestt, TaxaItem *mySourcet, int myType, int label){
    LWEdge *e=my_new(1, sizeof(LWEdge));
    e->source=mySourcet;
    e->dest=myDestt;
    e->type=myType;
    e->label=my_new(20, sizeof(char));
    sprintf(e->label, "%d", label);
    e->third=new_alist_l(16);
    e->inPath=0;
    e->path=NULL;
    return e;
}

void free_component(Component *cpt){
    return;

    int i, j;

    for(i=0;i<components->size;i++){
        Component *cur=components->elementData[i];
        if(cur==cpt) continue;
        Distance *minDist=cur->mindist;
        if(minDist->sc==cpt) return;
        if(minDist->dc==cpt) return;
    }

    if(cpt->mindist) free(cpt->mindist);
    if(cpt->compdist){
        //for(i=0;i<cpt->compdist->size;i++) free(cpt->compdist->elementData[i]);
        //free_alist_l(cpt->compdist);
    }
    if(cpt->taxa) free_alist_l(cpt->taxa);
    free(cpt);
}

//=================================

void calculate_connection_limit(){
    mylog("calculating connection limit...");

    int MaxNucleotides=seq_len;
    int conLimit=95;

    int index = 0;
    double parsval = 1;

    double percentValue = (conLimit * 0.01);

    while (parsval > percentValue) {
        index++;
        parsval = calcPars(index, MaxNucleotides - index, 2000);
    }

    maxParsimonyDistance=index-1;

    sprintf(loginfo, "maxParsimonyDistance=%d", maxParsimonyDistance);mylog(loginfo);
}

double calcPars(int j, int m, int it){
		double i;
		double product;
		double dq;
		double pq;
		double num1;
		double num2;
		double q;
		double int1;
		double int2;
		double t1;
		double tr;
		double tp;
		int b = 3;
		int r = 1;
		double u = 1.0;
		product = 1.0;
		dq = u / it;

		/*
		Product
		*/
		tr = (double)1.0 / (b * r);
		tp = (double)2.0 * m + 1.0;

		for (i = 0.00000001; i < j; i += 1.0) {

			/*
			Integrate
			*/
			int1 = 0.0;
			int2 = 0.0;

			for (q = 0.0000001; q < u; q += dq) {
				t1 = (double)1.0 - q * tr;
				num2 = pow(q, i) * pow(((double)1.0 - q), tp) * t1 * pow(t1 + 1.0 - q, i) * ((double)1.0 - 2 * q * t1) * dq;
				num1 = num2 * q;
				int1 += num1;
				int2 += num2;
			}

			pq = (double)1.0 - int1 / int2;
			product *= pq;
		}

		return (product);
}

void build_matrix(){
    mylog("building matrix...");

    if(1){
        split_span_for_threads(thread_num, 0, components->size);
        int i;
        for(i=0;i<thread_num;i++) thread_pool_add_worker(threadPool, build_matrix_one, i);
        thread_pool_invoke_all(threadPool);
    }else{
        int i, j, k;

        Distance *minDist=new_distance1();

        for(i=0;i<components->size;i++){
            Component *componentOne=components->elementData[i];
            TaxaItem *taxaItemOne=componentOne->taxa->elementData[0];
            ((Distance *)componentOne->mindist)->distance=INFINITY;
            minDist->distance=INFINITY;
            //--
            for(j=0;j<components->size;j++){
                Component *componentTwo=components->elementData[j];
                TaxaItem *taxaItemTwo=componentTwo->taxa->elementData[0];
                //--
                Distance *realdist=  new_distance3(componentOne->id, componentTwo->id, 0, componentOne, componentTwo);
                Distance *tmpcdist=  new_distance3(componentOne->id, componentTwo->id, 0, componentOne, componentOne);
                Distance *tmptdist=  new_distance3(taxaItemOne->id, taxaItemTwo->id, 0, componentOne, componentOne);
                Distance *metricdist=new_distance3(taxaItemOne->id, taxaItemTwo->id, 0, componentOne, componentOne);
                if(componentOne->id!=componentTwo->id){
                    for(k=0;k<taxaItemOne->l_characters;k++) {
                        char char1 = taxaItemOne->characters[k];
                        char char2 = taxaItemTwo->characters[k];
                        if (char1 == '?' || char2 == '?') continue;
                        if (char1 != char2) realdist->distance++;
                    }
                    tmpcdist->distance = realdist->distance;
                    tmptdist->distance = INFINITY;
                    metricdist->distance = INFINITY;
                    if(realdist->distance < minDist->distance){
                        clone_distance(minDist, realdist);
                    }
                }
                alist_l_add(componentOne->compdist, tmpcdist);
                alist_l_add(taxaItemOne->compdist, tmptdist);
                alist_l_add(taxaItemOne->realdist, realdist);
                alist_l_add(taxaItemOne->metricdist, metricdist);
            }
            clone_distance(componentOne->mindist, minDist);
            taxaItemOne->minRealDist = minDist->distance;
        }
        free(minDist);
    }
}

void build_matrix_one(int thread_index){
    int start=starts_thread[thread_index];
    int end=ends_thread[thread_index];

    int i, j, k;

    Distance *minDist=new_distance1();

    for(i=start;i<end;i++){
        Component *componentOne=components->elementData[i];
        TaxaItem *taxaItemOne=componentOne->taxa->elementData[0];
        ((Distance *)componentOne->mindist)->distance=INFINITY;
        minDist->distance=INFINITY;
        //--
        for(j=0;j<components->size;j++){
            Component *componentTwo=components->elementData[j];
            TaxaItem *taxaItemTwo=componentTwo->taxa->elementData[0];
            //--
            Distance *realdist=  new_distance3(componentOne->id, componentTwo->id, 0, componentOne, componentTwo);
            Distance *tmpcdist=  new_distance3(componentOne->id, componentTwo->id, 0, componentOne, componentOne);
            Distance *tmptdist=  new_distance3(taxaItemOne->id, taxaItemTwo->id, 0, componentOne, componentOne);
            Distance *metricdist=new_distance3(taxaItemOne->id, taxaItemTwo->id, 0, componentOne, componentOne);
            if(componentOne->id!=componentTwo->id){
                for(k=0;k<taxaItemOne->l_characters;k++) {
                    char char1 = taxaItemOne->characters[k];
                    char char2 = taxaItemTwo->characters[k];
                    if (char1 == '?' || char2 == '?') continue;
                    if (char1 != char2) realdist->distance++;
                }
                tmpcdist->distance = realdist->distance;
                tmptdist->distance = INFINITY;
                metricdist->distance = INFINITY;
                if(realdist->distance < minDist->distance){
                    clone_distance(minDist, realdist);
                }
            }
            alist_l_add(componentOne->compdist, tmpcdist);
            alist_l_add(taxaItemOne->compdist, tmptdist);
            alist_l_add(taxaItemOne->realdist, realdist);
            alist_l_add(taxaItemOne->metricdist, metricdist);
        }
        clone_distance(componentOne->mindist, minDist);
        taxaItemOne->minRealDist = minDist->distance;
    }
    free(minDist);
}

//=================================

int connectTaxa(TaxaItem *source){
    int madeNewConnection=0;

    int i, j;

    Distance *minDist=new_distance1();
    minDist->distance=INFINITY;

    for(i=0;i<source->realdist->size;i++){
        Distance *dist=source->realdist->elementData[i];
        if(dist->distance<minDist->distance && dist->destination!=source->id){
            clone_distance(minDist, dist);
        }
    }

    if(minDist->distance>maxParsimonyDistance){
        free(minDist);
        return madeNewConnection;
    }

    int is_combineComponents=0;

    MyIterator *iter=new_my_iterator(components);
    while(my_iterator_has_more(iter)){
        Component *curComponent=my_iterator_next(iter);
        for(i=0;i<curComponent->taxa->size;i++){
            TaxaItem *curTaxaItem=curComponent->taxa->elementData[i];
            if(curTaxaItem->isIntermediate) continue;
            Distance *dist=source->realdist->elementData[curTaxaItem->id];
            if(minDist->distance==dist->distance){
                is_combineComponents=1;
                Distance *already = source->compdist->elementData[dist->destination];
                if(already->distance <= minDist->distance) continue;
                int limit = source->minRealDist;
                if(limit<curTaxaItem->minRealDist) limit=curTaxaItem->minRealDist;
                Component *sourceParentComp=source->parentComponent;
                Component *curTaxaItemParentComp=curTaxaItem->parentComponent;
                AList_l *destCandidates=new_alist_l(curTaxaItemParentComp->taxa->size);
                for(j=0;j<curTaxaItemParentComp->taxa->size;j++) alist_l_add(destCandidates, curTaxaItemParentComp->taxa->elementData[j]);
                madeNewConnection|=bestMetric(sourceParentComp, curTaxaItemParentComp, source, curTaxaItem, limit, destCandidates);
                free_alist_l(destCandidates);
            }
        }
        if(is_combineComponents) {
            combineComponents(source->parentComponent, curComponent);
            free_component(curComponent);
            is_combineComponents=0;
        }
    }
    free(iter);
    free(minDist);

    return madeNewConnection;
}

int bestMetric(Component *sourcec, Component *destc, TaxaItem *sourcet, TaxaItem *destt, int limit, AList_l *destCandidates){
    int addedConnection=0;

    Distance *bestDist=new_distance2(0, 0, 0);
    int metric = -INFINITY;
    int firstTime = 1;

    int i, j, k;

    for(i=0;i<sourcec->taxa->size;i++){
        TaxaItem *curTaxaFromSrcsComp=sourcec->taxa->elementData[i];
        if(curTaxaFromSrcsComp!=sourcet && !curTaxaFromSrcsComp->isIntermediate) continue;
        for(j=0;j<destCandidates->size;j++){
            TaxaItem *curDestCandTaxa=destCandidates->elementData[j];
            if(curDestCandTaxa!=destt && !curDestCandTaxa->isIntermediate) continue;

            int total = ((Distance *)sourcet->realdist->elementData[destt->id])->distance;
            Distance *nearest = sourcet->compdist->elementData[curTaxaFromSrcsComp->id];
            Distance *other = destt->compdist->elementData[curDestCandTaxa->id];
            total -= nearest->distance + other->distance;
            if(total<1) continue;

            int srcCnt=sourcec->taxa->size;
            int dstCnt=0;
            AList_l *allvect=new_alist_l(srcCnt);
            for(k=0;k<srcCnt;k++) alist_l_add(allvect, sourcec->taxa->elementData[k]);
            if(sourcec->id!=destt->id){
                dstCnt=destc->taxa->size;
                for(k=0;k<dstCnt;k++) alist_l_add(allvect, destc->taxa->elementData[k]);
            }
            if(srcCnt>dstCnt) srcCnt*=10;
            else srcCnt=dstCnt*10;
            recalcDist(allvect, curTaxaFromSrcsComp, curDestCandTaxa, total);
            int score = evaluateMetric(allvect, srcCnt, total, limit);
            if(firstTime || (score>metric)) {
                firstTime=0;
                bestDist->source = curTaxaFromSrcsComp->id;
                bestDist->destination = curDestCandTaxa->id;
                bestDist->distance = total;
                metric = score;
            }
            free_alist_l(allvect);
        }
    }

    TaxaItem *bestDistSrcTaxa = alltaxa->elementData[bestDist->source];
    TaxaItem *bestDistDstTaxa = alltaxa->elementData[bestDist->destination];

    Distance *already = bestDistSrcTaxa->compdist->elementData[bestDistDstTaxa->id];

    if(already->distance > bestDist->distance){
        addedConnection = 1;
        connectComponents(sourcec, destc, bestDist, destCandidates);
        updateDistance(sourcec, destc, bestDist);
    }

    free(bestDist);

    return addedConnection;
}

int combineComponents(Component *componentOne, Component *componentTwo){
    if(componentOne==componentTwo || componentOne->id==componentTwo->id) return 0;

    int i;
    for(i=0;i<componentTwo->taxa->size;i++){
        TaxaItem *tempt1=componentTwo->taxa->elementData[i];
        alist_l_add(componentOne->taxa, tempt1);
        tempt1->parentComponent=componentOne;
    }

    remove_alist_l(components, componentTwo);
    recalcMinDistance(componentOne, componentTwo->id);

    return 1;
}

void print_taxa(TaxaItem *item){
    fprintf(stderr, "taxa_id=%d", item->id);
    if(item->parentComponent) fprintf(stderr, "\t%d", ((Component *)(item->parentComponent))->id);
    fprintf(stderr, "\n");

    int i;
    if(item->realdist){
        fprintf(stderr, "realdist:\n");
        for(i=0;i<item->realdist->size;i++) print_distance(item->realdist->elementData[i]);
    }
    if(item->compdist){
        fprintf(stderr, "compdist:\n");
        for(i=0;i<item->compdist->size;i++) print_distance(item->compdist->elementData[i]);
    }
    if(item->metricdist){
        fprintf(stderr, "metricdist:\n");
        for(i=0;i<item->metricdist->size;i++) print_distance(item->metricdist->elementData[i]);
    }
}

void print_distance(Distance *dis){
    fprintf(stderr, "distance: source=%d\tdest=%d\tdis=%d", dis->source, dis->destination, dis->distance);
    if(dis->sc) fprintf(stderr, "\tsc=%d", dis->sc->id);
    if(dis->dc) fprintf(stderr, "\tdc=%d", dis->dc->id);
    fprintf(stderr, "\n");
}

void print_component(Component *cpt){
    fprintf(stderr, "component_id=%d\ttaxa=%d\tmin_dis=%d\n", cpt->id, cpt->taxa->size, ((Distance *)(cpt->mindist))->distance);
    int i;
    for(i=0;i<cpt->taxa->size;i++) print_taxa(cpt->taxa->elementData[i]);
    if(cpt->compdist){
        for(i=0;i<cpt->compdist->size;i++) print_distance(cpt->compdist->elementData[i]);
    }
}

void remove_alist_l(AList_l *list, int64_t v){
    int i, size=list->size;

    int64_t *copy=my_new(size, sizeof(int64_t));
    memcpy(copy, list->elementData, size*sizeof(int64_t));

    list->size=0;
    for(i=0;i<size;i++){
        int64_t a=copy[i];
        if(a==v) continue;
        alist_l_add(list, a);
    }

    free(copy);
}

void recalcMinDistance(Component *collapse, int removeid){
    int i, j, k;

    for(i=0;i<components->size;i++){
        Component *tempc1=components->elementData[i];
        Distance *tmpDis=tempc1->mindist;
        if(tmpDis->dc->id==removeid) {
            tmpDis->dc=collapse;
        }
    }

    Distance *mindist = new_distance1();
    mindist->distance = INFINITY;

    for(i=0;i<components->size;i++){
        Component *tempc1=components->elementData[i];
        if(tempc1==collapse) continue;
        for(j=0;j<tempc1->taxa->size;j++){
            TaxaItem *tempt1=tempc1->taxa->elementData[j];
            if(tempt1->isIntermediate) continue;
            for(k=0;k<collapse->taxa->size;k++){
                TaxaItem *tempt2=collapse->taxa->elementData[k];
                if(tempt2->isIntermediate) continue;
                Distance *nextdist = tempt2->realdist->elementData[tempt1->id];
                if (nextdist->distance < mindist->distance){
                    clone_distance(mindist, nextdist);
                    mindist->sc = collapse;
                    mindist->dc = tempc1;
                }
            }
        }
    }

    if(collapse->mindist) free(collapse->mindist);
    collapse->mindist=mindist;
}

TaxaItem *recalcDist_sourcebordert;
TaxaItem *recalcDist_destbordert;
int recalcDist_length;

void recalcDist(AList_l *allvect, TaxaItem *sourcebordert, TaxaItem *destbordert, int length){
    Distance *metricDist = sourcebordert->metricdist->elementData[destbordert->id];
    metricDist->distance = length;

    int i, j;

    allvect_for_multi_thread=allvect;
    recalcDist_sourcebordert=sourcebordert;
    recalcDist_destbordert=destbordert;
    recalcDist_length=length;
    for(i=0;i<allvect_for_multi_thread->size;i++){
        thread_pool_add_worker(threadPool, recalcDist_one, i);
        //recalcDist_one(i);
    }
    thread_pool_invoke_all(threadPool);
}

void recalcDist_one(int64_t para){
    int thread_index=(int)(para>>32);
    int j, i=(int)(para&0xFFFFFFFFL);

    TaxaItem *sourcebordert=recalcDist_sourcebordert;
    TaxaItem *destbordert=recalcDist_destbordert;
    int length=recalcDist_length;

    {
        TaxaItem *sourceitem = allvect_for_multi_thread->elementData[i];
        for(j=0;j<allvect_for_multi_thread->size;j++){
            TaxaItem *destitem = allvect_for_multi_thread->elementData[j];
            //--
            Distance *dist = sourceitem->metricdist->elementData[destitem->id];
            Distance *compdist = sourceitem->compdist->elementData[destitem->id];
            int newdistance = ((Distance *)sourceitem->compdist->elementData[sourcebordert->id])->distance;

            newdistance += length;
            newdistance += ((Distance *)destbordert->compdist->elementData[destitem->id])->distance;

            if (newdistance < compdist->distance) dist->distance = newdistance;
            else dist->distance = compdist->distance;
        }
    }
}

int evaluateMetric_TOO_SMALL_SCORE;
int evaluateMetric_limit;

int evaluateMetric(AList_l *allvect, int TOO_SMALL_SCORE, int numIntermediates, int limit){
    int score = -(numIntermediates - 1);

    int i, j;

    AList_l *list=new_alist_l(16);
    for(i=0;i<allvect->size;i++){
        TaxaItem *item=allvect->elementData[i];
        if(!item->isIntermediate) alist_l_add(list, item);
    }

    allvect_for_multi_thread=list;
    set_evaluateMetric_is_inf(0);
    memset(evaluateMetric_results, 0, thread_num*sizeof(int));
    evaluateMetric_TOO_SMALL_SCORE=TOO_SMALL_SCORE;
    evaluateMetric_limit=limit;
    for(i=0;i<allvect_for_multi_thread->size;i++){
        thread_pool_add_worker(threadPool, evaluateMetric_one, i);
    }
    thread_pool_invoke_all(threadPool);

    if(get_evaluateMetric_is_inf()){
        score=-INFINITY;
    }else{
        for(i=0;i<thread_num;i++) score+=evaluateMetric_results[i];
    }

    free(list);

    return score;
}

void evaluateMetric_one(int64_t para){
    int thread_index=(int)(para>>32);
    int j, i=(int)(para&0xFFFFFFFFL);

    int TOO_SMALL_SCORE=evaluateMetric_TOO_SMALL_SCORE;
    int limit=evaluateMetric_limit;

    int score=0;
    {
        if(get_evaluateMetric_is_inf()){score=-INFINITY;goto last;}
        TaxaItem *sourceitem = allvect_for_multi_thread->elementData[i];
        for(j=0;j<allvect_for_multi_thread->size;j++){
            TaxaItem *destitem = allvect_for_multi_thread->elementData[j];
            if(sourceitem->id==destitem->id) continue;
            //--
            Distance *realdist = sourceitem->realdist->elementData[destitem->id];
            Distance *metricdist = sourceitem->metricdist->elementData[destitem->id];

            if((realdist->distance > (2 * maxParsimonyDistance)) && (metricdist->distance > (2 * maxParsimonyDistance))) continue;
            if(realdist->distance == metricdist->distance) {
                score += CORRECT_SCORE;
            }else{
                if(realdist->distance > metricdist->distance) {
                    if((sourceitem->parentComponent != destitem->parentComponent) && (metricdist->distance <= limit)){
                        score = -INFINITY;
                        set_evaluateMetric_is_inf(1);
                        goto last;
                    }else{
                        double exp_factor=realdist->distance - metricdist->distance;
                        score -= TOO_SMALL_SCORE * exp_factor;
                    }
                }else{
                    score -= TOO_BIG_SCORE;
                }
            }
            if(metricdist->distance < sourceitem->minRealDist){
                score = -INFINITY;
                set_evaluateMetric_is_inf(1);
                goto last;
            }
        }
    }

    last:
    return evaluateMetric_results[thread_index]+=score;
}

void set_evaluateMetric_is_inf(int is_inf){
    pthread_mutex_lock(locker);
    evaluateMetric_is_inf=is_inf;
    pthread_mutex_unlock(locker);
}

int get_evaluateMetric_is_inf(){
    int is_inf;
    pthread_mutex_lock(locker);
    is_inf=evaluateMetric_is_inf;
    pthread_mutex_unlock(locker);
    return is_inf;
}

void connectComponents(Component *sourcec, Component *destc, Distance *newdist, AList_l *candidates){
    TaxaItem *sourcet = alltaxa->elementData[newdist->source];
    TaxaItem *destt = alltaxa->elementData[newdist->destination];

    graphExists = 1;

    if(newdist->distance == 1){
        alist_l_add(sourcet->nbor, destt);
        alist_l_add(destt->nbor, sourcet);
        addLWEdge(sourcet, destt);
    }else{
        TaxaItem *nextt = sourcet;
        int intermediates = newdist->distance;
        char *consensus = NULL;
        int l_consensus=0;

        while(intermediates>1){
            intermediates--;

            sprintf(loginfo, "IN%d", internalNodeNumber);
            TaxaItem *newt = new_taxaitem1(loginfo, 0, consensus, l_consensus);
            newt->isIntermediate = 1;

            internalNodeNumber++;

            alist_l_add(alltaxa, newt);
            alist_l_add(candidates, newt);
            newt->id=alltaxa->size-1;

            alist_l_add(sourcec->taxa, newt);

            alist_l_add(nextt->nbor, newt);
            alist_l_add(newt->nbor, nextt);
            //fprintf(stderr, "1\t%d\t%d\n", nextt->id, newt->id);
            addLWEdge(nextt, newt);

            addIntermediate(newt, nextt, destt, intermediates);

            nextt = newt;
            if(intermediates == 1) {
                alist_l_add(destt->nbor, newt);
                alist_l_add(newt->nbor, destt);
                //fprintf(stderr, "2\t%d\t%d\n", destt->id, newt->id);
                addLWEdge(newt, destt);
                break;
            }
        }
    }
}

void updateDistance(Component *sourcec, Component *destc, Distance *newdist){
    int i, j;
    for(i=0;i<sourcec->taxa->size;i++){
        TaxaItem *curSrcTaxa=sourcec->taxa->elementData[i];
        Distance *curSrcDist=curSrcTaxa->compdist->elementData[newdist->source];
        for(j=0;j<destc->taxa->size;j++){
            TaxaItem *curDestTaxa=destc->taxa->elementData[j];
            Distance *curDestDist=curDestTaxa->compdist->elementData[newdist->destination];
            if((curSrcDist->distance==INFINITY) || (curDestDist->distance==INFINITY)) continue;
            else{
                Distance *curDist=curSrcTaxa->compdist->elementData[curDestTaxa->id];
                int newDistance=curDestDist->distance+curSrcDist->distance+newdist->distance;
                if(newDistance < curDist->distance){
                    curDist->distance = newDistance;
                    curDist=curDestTaxa->compdist->elementData[curSrcTaxa->id];
                    curDist->distance = newDistance;
                }
            }
        }
    }
}

void addLWEdge(TaxaItem *sourcet, TaxaItem *destt){
    if(destt->isIntermediate ^ sourcet->isIntermediate){
        alist_l_add(edges, new_LWEdge(destt, sourcet, LWEdge_HAPINT, edges->size));
        numHapInt++;
    }else if (!destt->isIntermediate && !sourcet->isIntermediate){
        alist_l_add(edges, new_LWEdge(destt, sourcet, LWEdge_HAPHAP, edges->size));
        numHapHap++;
    }else{
        alist_l_add(edges, new_LWEdge(destt, sourcet, LWEdge_INTINT, edges->size));
    }
}

void addIntermediate(TaxaItem *newtaxa, TaxaItem *before, TaxaItem *after, int after_distance){
    int i;
    for(i=0;i<alltaxa->size;i++){
        TaxaItem *curtaxa=alltaxa->elementData[i];
        Distance *newdistance;
        Distance *fromnewdistance;
        Distance *dummy;

        if(curtaxa==newtaxa) {
            fromnewdistance = new_distance2(newtaxa->id, curtaxa->id, 0);
            alist_l_add(newtaxa->compdist, fromnewdistance);
            dummy = new_distance2(newtaxa->id, curtaxa->id, 0);
            alist_l_add(newtaxa->metricdist, dummy);
            continue;
        }

        Distance *beforetaxa = curtaxa->compdist->elementData[before->id];
        Distance *aftertaxa = curtaxa->compdist->elementData[after->id];

        if((beforetaxa->distance != INFINITY) && (beforetaxa->distance < aftertaxa->distance)){
            newdistance = new_distance2(curtaxa->id, newtaxa->id, beforetaxa->distance+1);
            fromnewdistance = new_distance2(newtaxa->id, curtaxa->id, beforetaxa->distance+1);
        }else if(aftertaxa->distance != INFINITY){
            newdistance = new_distance2(curtaxa->id, newtaxa->id, aftertaxa->distance+after_distance);
            fromnewdistance = new_distance2(newtaxa->id, curtaxa->id, aftertaxa->distance+after_distance);
        }else{
            newdistance = new_distance2(curtaxa->id, newtaxa->id, INFINITY);
            fromnewdistance = new_distance2(newtaxa->id, curtaxa->id, INFINITY);
        }

        alist_l_add(curtaxa->compdist, newdistance);
        alist_l_add(newtaxa->compdist, fromnewdistance);
        dummy = new_distance2(curtaxa->id, newtaxa->id, INFINITY);
        alist_l_add(curtaxa->metricdist, dummy);
        dummy = new_distance2(newtaxa->id, curtaxa->id, INFINITY);
        alist_l_add(newtaxa->metricdist, dummy);
    }
}

int connectComponents2(){
    Distance *mindist = new_distance1();
    mindist->distance = INFINITY;

    int i, j;

    for(i=0;i<components->size;i++){
        Component *tempc1=components->elementData[i];
        Distance *tmpDis=tempc1->mindist;
        if(tmpDis->distance < mindist->distance) clone_distance(mindist, tmpDis);
    }

    if(mindist->distance>maxParsimonyDistance){free(mindist);return 0;}

    Component *tempc1 = mindist->sc;
    Component *tempc2 = mindist->dc;

    AList_l *destCandidates=new_alist_l(mindist->dc->taxa->size);
    for(i=0;i<mindist->dc->taxa->size;i++) alist_l_add(destCandidates, mindist->dc->taxa->elementData[i]);

    MyIterator *iter1=new_my_iterator(tempc1->taxa);
    while(my_iterator_has_more(iter1)){
        TaxaItem *tempt1=my_iterator_next(iter1);
        if(tempt1->isIntermediate) continue;
        MyIterator *iter2=new_my_iterator(tempc2->taxa);
        while(my_iterator_has_more(iter2)){
            TaxaItem *tempt2=my_iterator_next(iter2);
            if(tempt2->isIntermediate) continue;
            Distance *thisdist = tempt1->realdist->elementData[tempt2->id];
            if(thisdist->distance == mindist->distance) {
                Distance *curdist = tempt1->compdist->elementData[tempt2->id];
                if(curdist->distance <= mindist->distance) continue;
                bestMetric(tempc1, tempc2, tempt1, tempt2, mindist->distance, destCandidates);
            }
        }
        free(iter2);
    }
    free(iter1);
    free_alist_l(destCandidates);
    free(mindist);

    combineComponents(tempc1, tempc2);
    free_component(tempc2);

    return 1;
}

//=================================

void calculateOutgroupWeights(){
    int numNetworks = components->size + 1;
    double *total_weight=my_new(numNetworks, sizeof(double));
    int networkNum = 0;
    int ntaxa = 0;

    int i, j, k;

    for(i=0;i<components->size;i++){
        Component *tempc=components->elementData[i];
        networkNum++;
        for(j=0;j<tempc->taxa->size;j++){
            TaxaItem *tempt=tempc->taxa->elementData[j];
            if(tempt->isIntermediate)  continue;
            ntaxa++;
            int sum_nbor_dup = 0;
            int anynbor = 0;
            for(k=0;k<tempt->nbor->size;k++){
                TaxaItem *tempn=tempt->nbor->elementData[k];
                if(!tempn->isIntermediate) sum_nbor_dup += (tempn->numduplicates+1);
                anynbor++;
            }
            if(anynbor==1) tempt->isTip=1;
            else tempt->isTip=0;
            if(tempt->isTip) tempt->weight = (tempt->numduplicates+1)/2.0;
            else tempt->weight = tempt->numduplicates+1+sum_nbor_dup;
            total_weight[networkNum] += tempt->weight;
        }
    }

    maxWtaxa = my_new(numNetworks, sizeof(char *));
    double *maxW = my_new(numNetworks, sizeof(double));
    networkNum = 0;
    maxW[networkNum] = 0.0;

    for(i=0;i<components->size;i++){
        Component *tempcomp=components->elementData[i];
        networkNum++;
        maxW[networkNum] = 0.0;
        for(j=0;j<tempcomp->taxa->size;j++){
            TaxaItem *temptaxa=tempcomp->taxa->elementData[j];
            if(temptaxa->isIntermediate) continue;
            temptaxa->oweight = temptaxa->weight / total_weight[networkNum];
            if(temptaxa->oweight > maxW[networkNum]){
                maxW[networkNum] = temptaxa->oweight;
                maxWtaxa[networkNum] = temptaxa->name;
            }
        }
    }

    free(maxW);
    free(total_weight);
}

void get_LWEdge_diff(LWEdge *edge, char *res_c1, char *res_c2){
    TaxaItem *source=edge->source;
    TaxaItem *target=edge->dest;

    int i;
    for(i=0;i<source->l_characters;i++){
        char c1=source->characters[i];
        char c2=target->characters[i];
        if(c1!=c2){
            *(res_c1)=c1;
            *(res_c2)=c2;
            break;
        }
    }
}

void printGraph(){
    int i, j, k, len;
    char c1, c2;

    sprintf(loginfo, "%s.gml", outFile);
    GzStream *out=gz_stream_open(loginfo, "w");

    gz_write(out, "graph [\n", 8);
    gz_write(out, "   directed 0\n", 14);

    for(i=0;i<components->size;i++){
        Component *cpt=components->elementData[i];
        for(j=0;j<cpt->taxa->size;j++){
            TaxaItem *item=cpt->taxa->elementData[j];
            //len=sprintf(loginfo, "vertex:\t%d\t%s\t%d\t%f\n", item->id, item->name, item->isIntermediate, item->oweight);
            //gz_write(out2, loginfo, len);
            gz_write(out, "   node [\n", 10);
            gz_write(out, "      data [\n", 13);
            gz_write(out, loginfo, sprintf(loginfo, "         Frequency \"frequency=%.1f\n", item->isIntermediate ? 1.0:(double)(item->dupnames->size+1)));
            gz_write(out, loginfo, sprintf(loginfo, "%s\n", item->name));
            for(k=0;k<item->dupnames->size;k++) gz_write(out, loginfo, sprintf(loginfo, "%s\n", item->dupnames->elementData[k]));
            gz_write(out, "\"\n      ]\n", 10);
            gz_write(out, loginfo, sprintf(loginfo, "      id %d\n   ]\n", item->id));
        }
    }

    for(i=0;i<edges->size;i++){
        LWEdge *edge=edges->elementData[i];
        //len=sprintf(loginfo, "edge:\t%d\t%d\t1\n", edge->source->id, edge->dest->id);
        //gz_write(out2, loginfo, len);
        //get_LWEdge_diff(edge, &c1, &c2);
        gz_write(out, "   edge [\n", 10);
        gz_write(out, "      data [\n", 13);
        gz_write(out, loginfo, sprintf(loginfo, "         Changes \"1\"\n"));
        gz_write(out, "      ]\n", 8);
        gz_write(out, loginfo, sprintf(loginfo, "      source %d\n", edge->source->id));
        gz_write(out, loginfo, sprintf(loginfo, "      target %d\n", edge->dest->id));
        gz_write(out, "   ]\n", 5);
    }

    gz_write(out, "]\n", 2);

    gz_stream_destory(out);

    sprintf(loginfo, "%s.json", outFile);
    out=gz_stream_open(loginfo, "w");

    gz_write(out, "{\n", 2);

    gz_write(out, "  \"nodes\": [\n", 13);
    for(i=0;i<components->size;i++){
        Component *cpt=components->elementData[i];
        for(j=0;j<cpt->taxa->size;j++){
            TaxaItem *item=cpt->taxa->elementData[j];
            SBuilder *sb=new_s_builder(16);
            s_builder_add_str(sb, item->name);
            for(k=0;k<item->dupnames->size;k++){
                s_builder_add_char(sb, ';');
                s_builder_add_str(sb, item->dupnames->elementData[k]);
            }
            double freq=item->isIntermediate ? 1.0:(double)(item->dupnames->size+1);
            gz_write(out, loginfo, sprintf(loginfo, "    {\"id\": %d, \"frequency\": %f, \"title1\": \"%s\", \"title2\": \"%s\"}", item->id, freq, item->name, sb->str));
            if(i<(components->size-1) || j<(cpt->taxa->size-1)) gz_write_char(out, ',');
            gz_write_char(out, '\n');
            free_s_builder(sb);
        }
    }
    gz_write(out, "  ],\n", 5);

    gz_write(out, "  \"edges\": [\n", 13);
    for(i=0;i<edges->size;i++){
        LWEdge *edge=edges->elementData[i];
        gz_write(out, loginfo, sprintf(loginfo, "    {\"source\": %d, \"target\": %d, \"change\": 1}", edge->source->id, edge->dest->id));
        if(i<(edges->size-1)) gz_write_char(out, ',');
        gz_write_char(out, '\n');
    }
    gz_write(out, "  ]\n", 4);

    gz_write(out, "}\n", 2);

    gz_stream_destory(out);
}

//=================================

void get_TCSEdge_diff(TCSEdge *edge, char *res_c1, char *res_c2){
    TCSNode *source=edge->from;
    TCSNode *target=edge->to;

    int len=strlen(source->seq);

    int i;
    for(i=0;i<len;i++){
        char c1=source->seq[i];
        char c2=target->seq[i];
        if(c1!=c2){
            *(res_c1)=c1;
            *(res_c2)=c2;
            break;
        }
    }
}

void output_graph2(){
    mylog("constructing new graph...");

    int i, j, k, len;
    char c1, c2;

    int l_nodes=10000;
    int64_t *all_nodes=my_new(l_nodes, sizeof(int64_t));
    AList_l *all_edges=new_alist_l(16);

    for(i=0;i<components->size;i++){
        Component *cpt=components->elementData[i];
        for(j=0;j<cpt->taxa->size;j++){
            TaxaItem *item=cpt->taxa->elementData[j];
            //--
            TCSNode *node=my_new(1, sizeof(TCSNode));
            node->id=item->id;
            node->name=item->name;
            node->seq=item->characters;
            node->freq=item->numduplicates+1;
            node->is_intermediate=item->isIntermediate;
            node->dupnames=item->dupnames;
            node->in_edges=new_llist_l();
            node->out_edges=new_llist_l();
            node->is_output=1;
            //--
            if(item->id>=l_nodes){
                int span=item->id+100-l_nodes;
                all_nodes=my_renew(all_nodes, (l_nodes+span)*sizeof(int64_t));
                memset(all_nodes+l_nodes, 0, span*sizeof(int64_t));
                l_nodes+=span;
            }
            all_nodes[item->id]=node;
        }
    }

    Hash_li *hp=new_hash_li1(edges->size*2);
    for(i=0;i<edges->size;i++){
        LWEdge *edge=edges->elementData[i];
        TCSNode *from=all_nodes[edge->source->id];
        TCSNode *to=all_nodes[edge->dest->id];
        //--
        int64_t key1=from->id;
        int64_t key2=to->id;
        int64_t key3=key1<<32|key2;
        int64_t key4=key2<<32|key1;
        if(hash_li_get(hp, key3)) continue;
        if(hash_li_get(hp, key4)) continue;
        hash_li_put(hp, key3, 1);
        hash_li_put(hp, key4, 1);
        //--
        TCSEdge *e=my_new(1, sizeof(TCSEdge));
        e->from=from;
        e->to=to;
        e->weight=1;
        e->is_output=1;
        //--
        llist_l_add(from->out_edges, e);
        llist_l_add(to->in_edges, e);
        alist_l_add(all_edges, e);
    }
    free_hash_li(hp);

    mylog("merging graph intermediate nodes...");
    int is_merge;
    do{
        is_merge=0;
        for(i=0;i<l_nodes;i++){
            TCSNode *node=all_nodes[i];
            if(!node) continue;
            if(!node->is_intermediate) continue;
            if(!node->is_output) continue;
            //--
            if(node->in_edges->size==1 && node->out_edges->size==1){
                is_merge=1;
                node->is_output=0;
                //--
                TCSEdge *in=node->in_edges->first->value;
                TCSEdge *out=node->out_edges->first->value;
                TCSNode *from=in->from;
                TCSNode *to=out->to;
                //--
                TCSEdge *e=my_new(1, sizeof(TCSEdge));
                e->from=from;
                e->to=to;
                e->weight=in->weight+out->weight;
                e->is_output=1;
                //--
                llist_l_add(from->out_edges, e);
                llist_l_add(to->in_edges, e);
                alist_l_add(all_edges, e);
                //--
                llist_l_remove2(from->out_edges, in);
                llist_l_remove2(to->in_edges, out);
                in->is_output=0;
                out->is_output=0;
            }
        }
    }while(is_merge);

    mylog("outputting graph...");

    sprintf(loginfo, "%s.gml", outFile);
    GzStream *out=gz_stream_open(loginfo, "w");

    gz_write(out, "graph [\n", 8);
    gz_write(out, "   directed 0\n", 14);

    for(i=0;i<l_nodes;i++){
        TCSNode *node=all_nodes[i];
        if(!node) continue;
        if(!node->is_output) continue;
        //len=sprintf(loginfo, "vertex:\t%d\t%s\t%d\t%d\n", node->id, node->name, node->is_intermediate, node->freq);
        //gz_write(out2, loginfo, len);
        gz_write(out, "   node [\n", 10);
        gz_write(out, "      data [\n", 13);
        gz_write(out, loginfo, sprintf(loginfo, "         Frequency \"frequency=%.1f\n", node->is_intermediate ? 1.0:(double)(node->freq)));
        gz_write(out, loginfo, sprintf(loginfo, "%s\n", node->name));
        for(k=0;k<node->dupnames->size;k++) gz_write(out, loginfo, sprintf(loginfo, "%s\n", node->dupnames->elementData[k]));
        gz_write(out, "\"\n      ]\n", 10);
        gz_write(out, loginfo, sprintf(loginfo, "      id %d\n   ]\n", node->id));
    }

    for(i=0;i<all_edges->size;i++){
        TCSEdge *edge=all_edges->elementData[i];
        if(!edge->is_output) continue;
        //len=sprintf(loginfo, "edge:\t%d\t%d\t%d\n", edge->from->id, edge->to->id, edge->weight);
        //gz_write(out2, loginfo, len);
        //get_TCSEdge_diff(edge, c1, c2);
        gz_write(out, "   edge [\n", 10);
        gz_write(out, "      data [\n", 13);
        gz_write(out, loginfo, sprintf(loginfo, "         Changes \"1\"\n"));
        gz_write(out, "      ]\n", 8);
        gz_write(out, loginfo, sprintf(loginfo, "      source %d\n", edge->from->id));
        gz_write(out, loginfo, sprintf(loginfo, "      target %d\n", edge->to->id));
        gz_write(out, "   ]\n", 5);
    }

    gz_write(out, "]\n", 2);

    gz_stream_destory(out);

    sprintf(loginfo, "%s.json", outFile);
    out=gz_stream_open(loginfo, "w");

    gz_write(out, "{\n", 2);

    gz_write(out, "  \"nodes\": [\n", 13);
    int last_node_index=0;
    for(i=0;i<l_nodes;i++){
        TCSNode *node=all_nodes[i];
        if(!node) continue;
        if(!node->is_output) continue;
        last_node_index=i;
    }
    for(i=0;i<l_nodes;i++){
        TCSNode *node=all_nodes[i];
        if(!node) continue;
        if(!node->is_output) continue;
        SBuilder *sb=new_s_builder(16);
        s_builder_add_str(sb, node->name);
        for(k=0;k<node->dupnames->size;k++){
            s_builder_add_char(sb, ';');
            s_builder_add_str(sb, node->dupnames->elementData[k]);
        }
        double freq=node->is_intermediate ? 1.0:(double)(node->freq);
        gz_write(out, loginfo, sprintf(loginfo, "    {\"id\": %d, \"frequency\": %f, \"title1\": \"%s\", \"title2\": \"%s\"}", node->id, freq, node->name, sb->str));
        if(i<last_node_index) gz_write_char(out, ',');
        gz_write_char(out, '\n');
        free_s_builder(sb);
    }
    gz_write(out, "  ],\n", 5);

    gz_write(out, "  \"edges\": [\n", 13);
    for(i=0;i<all_edges->size;i++){
        TCSEdge *edge=all_edges->elementData[i];
        if(!edge->is_output) continue;
        gz_write(out, loginfo, sprintf(loginfo, "    {\"source\": %d, \"target\": %d, \"change\": 1}", edge->from->id, edge->to->id));
        if(i<(all_edges->size-1)) gz_write_char(out, ',');
        gz_write_char(out, '\n');
    }
    gz_write(out, "  ]\n", 4);

    gz_write(out, "}\n", 2);

    gz_stream_destory(out);
}

//=================================

