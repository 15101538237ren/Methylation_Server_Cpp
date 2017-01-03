//
//  main.cpp
//  Methylation_Server
//
//  Created by 任红雷 on 2016/12/29.
//  Copyright (c) 2016 ___renhonglei___. All rights reserved.
//

#include<iostream>
#include<cstdio>
#include<cmath>
#include<sstream>
#include<string>
#include<vector>
#include <algorithm>
#include <functional>
#include<cstdlib>
#include<ctime>
#include<time.h>
#include<regex>
#include <numeric>
#include <functional>
#include <fstream>
#include <iomanip>
#include <map>
#include <random>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#define pb push_back
using namespace std;
const float EPSINON = 1e-6;

double beta=0.195,gamma_param=0.04,power=1.0,alpha=12.0,c_off=0.13;
double pi = acos(-1);
int propencity_list = 9;
char right_status_of_reaction[] = {'U','H','M','H','H','H','U','H','M'};
char right_status_hash[] = {'H','M','H','U','M','M','H','U','H'};
int exist_or_make_dir(string dir_path)
{
    if (access(dir_path.c_str(), 0) == -1)
    {
        cout<<dir_path<<" is not existing"<<endl;
        cout<<"now make it"<<endl;
        int flag=mkdir(dir_path.c_str(), 0777);
        if (flag == 0)
        {
            cout<<"make successfully"<<endl;
            return 0;
        } else {
            cout<<"make errorly"<<endl;
            return 1;
        }
    }
    return 0;
}

int exist_or_delete_dir(string dir)
{
    if (access(dir.c_str(), 0) == 0)
    {
        cout<<dir<<" exists"<<endl;
        cout<<"now delete it"<<endl;
        int flag=rmdir(dir.c_str());
        if (flag == 0)
        {
            cout<<"delete it successfully"<<endl;
            return 0;
        } else {
            cout<<"delete it errorly"<<endl;
            return 1;
        }
    }
    return 0;
}

vector<int> index_random(int sample_num, int range){
    vector<int> index,return_index;
    srand(time(0));
    for(int i=0; i<range; i++){
        index.pb(0);
    }
    for(int i=1; i<=sample_num; i++){
        int idx = rand()%range;
        if(index[idx]==0){
            index[idx] = 1;
        }else{
            while(1){
                idx = rand()%range;
                if(index[idx]==0){
                    index[idx] = 1;
                    break;
                }
            }
        }
    }
    for(int i=0; i<range; i++){
        if(index[i]==1){
            return_index.pb(i);
        }
    }
    return return_index;
}

double phi_func(int distance){
    
    double fx = (1.0 / pow((beta + gamma_param * distance),power));
    double sx = (alpha / (30.0 * sqrt(2 * pi))) * exp(-(distance - 160.0) * (distance - 160.0) / (2 * 30.0 * 30.0));
    double rtn_phi = fx + sx + c_off;
    return rtn_phi;
    
}

vector<double> calc_propensity_list(double phi_d, vector<double> propensity_list, double pij, char xj_status){
    double U_plus=propensity_list[0];
    double H_plus=propensity_list[1];
    double M_minus=propensity_list[2];
    double H_minus=propensity_list[3];
    double H_p_H=propensity_list[4];
    double H_p_M=propensity_list[5];
    double U_p_M=propensity_list[6];
    double H_m_U=propensity_list[7];
    double M_m_U=propensity_list[8];
    
    double u_i_plus=U_plus+pij*double(xj_status=='M')*phi_d*(U_p_M-U_plus);
    double h_i_plus=H_plus+pij*double(xj_status=='M')*phi_d*(H_p_M-H_plus) + pij*double(xj_status=='H')*phi_d*(H_p_H-H_plus);
    double m_i_minus=M_minus+pij*double(xj_status=='U')*phi_d*(M_m_U-M_minus);
    double h_i_minus=H_minus+pij*double(xj_status=='U')*phi_d*(H_m_U-H_minus);
    vector<double> rtn_list = {u_i_plus,h_i_plus,m_i_minus,h_i_minus};
    return rtn_list;
    
}

int select_reaction(vector<double> propensity_tmp, int num_of_reactions, double sum_propensity, double random_number){
    int reaction = -1;
    double tmp_sum_propencity = 0.0;
    random_number = random_number * sum_propensity;
    for(int i = 0; i<num_of_reactions; i++){
        tmp_sum_propencity = tmp_sum_propencity + propensity_tmp[i];
        if(random_number < tmp_sum_propencity){
            reaction = i;
            break;
        }
    }
    return reaction;
}

void simulate(int round_no,int round_start,int out_target_start_gen,int out_target_end_gen,int generation, int time_step, string init_cell, string detail_file_dir, string ratio_file_dir, vector<double> propensity_list, vector<int> index_pos_list, int nearby = -1, int max_cells = 1000, bool real_nearby = false){
    
    
    vector<string> cell_collection;
    cell_collection.pb(init_cell);
    srand(time(0));
    for(int i=1; i<=generation; i++){
        clock_t gen_start_time=clock();
        time_t raw_time_start;
        struct tm * start_time_info;
        time ( &raw_time_start );
        start_time_info = localtime ( &raw_time_start );
        printf ( "round %d, generataion %d start time is: %s", round_no, i,asctime (start_time_info) );
        
        
        vector<string> cells_wait_to_add;
        unsigned long cc_size = cell_collection.size();
        if(cc_size > max_cells/2){
            vector<string>::iterator it=cell_collection.begin();
            vector<int> index_vec = index_random(max_cells/2, max_cells);
            int pos=0;
            for(int j=0;it!=cell_collection.end();j++){
                if((pos>=index_vec.size()) or (j != index_vec[pos])){
                    cell_collection.erase(it);
                }else{
                    pos++;
                    it++;
                }
            }
        }
        
        vector<vector<int> > M_count_statistics(time_step);
        vector<vector<int> > H_count_statistics(time_step);
        vector<vector<int> > U_count_statistics(time_step);
        
        vector<vector<string> > out_detail_seq_arr(cell_collection.size());
        for(int idx=0; idx<cell_collection.size(); idx++){
            vector<string> vect_tmp;
            for(int j=0; j<time_step; j++){
                clock_t ts_start_time=clock();
                int cell_len = cell_collection[idx].length();
                
                for(int k=0; k<cell_len; k++){
                    int target_reaction_CpG_site = 1 + rand()%(cell_len-2);  //1~cell_len-2
                    int col_CpG_site_index = target_reaction_CpG_site + ((rand()%2)?1:-1);
                    char status_col_site = cell_collection[idx][col_CpG_site_index];
                    
                    int pos_target = index_pos_list[target_reaction_CpG_site];
                    int col_site_pos = index_pos_list[col_CpG_site_index];
                    int distance = int(fabs(pos_target - col_site_pos));
                    
                    if(real_nearby == true and distance > nearby){
                        continue;
                    }
                    if(distance > 997){
                        continue;
                    }
                    double phi_d = phi_func(distance);
                    double pij = 1.0;
                    vector<double> propensity_tmp = calc_propensity_list(phi_d, propensity_list, pij, status_col_site);
                    double sum_propensity = 0.0;
                    int num_of_reactions = propensity_tmp.size();
                    for(int i=0; i<num_of_reactions; i++){
                        sum_propensity += propensity_tmp[i];
                    }
                    
                    double random_num = (double)(rand()%RAND_MAX)/RAND_MAX;
                    int reaction_id = select_reaction(propensity_tmp, num_of_reactions, sum_propensity, random_num);
                    char status_of_target_site = cell_collection[idx][target_reaction_CpG_site];
                    if(right_status_of_reaction[reaction_id] != status_of_target_site){
                        continue;
                    }
                    
                    cell_collection[idx][target_reaction_CpG_site] = right_status_hash[reaction_id];
                }
                int m_count=0,h_count=0,u_count=0;
                for (int it=0;it < cell_len;it++)
                {
                    switch (cell_collection[idx][it]) {
                        case 'M':
                            m_count++;
                            break;
                        case 'H':
                            h_count++;
                            break;
                        case 'U':
                            u_count++;
                            break;
                        default:
                            break;
                    }
                }
                M_count_statistics[j].pb(m_count);
                H_count_statistics[j].pb(h_count);
                U_count_statistics[j].pb(u_count);
                
                out_detail_seq_arr[idx].pb(cell_collection[idx]);
                clock_t ts_end_time=clock();
                //cout<< "gen:"<< i << ", time step:" << j << ": "<<static_cast<double>(ts_end_time-ts_start_time)/CLOCKS_PER_SEC*1000<<"ms"<<endl;//输出运行时间
            }
            
            string cell_1, cell_2;
            for(int j=0; j<cell_collection[idx].length(); j++){
                char ch = cell_collection[idx][j];
                if(ch == 'M'){
                    cell_1.append(1, 'H');
                    cell_2.append(1, 'H');
                }else if(ch == 'U'){
                    cell_1.append(1, 'U');
                    cell_2.append(1, 'U');
                }else if(ch == 'H'){
                    cell_1.append(1, 'H');
                    cell_2.append(1, 'U');
                }
            }
            cells_wait_to_add.pb(cell_1);
            cells_wait_to_add.pb(cell_2);
        }
        
        vector<double> m_means_ratio;
        vector<double> h_means_ratio;
        vector<double> u_means_ratio;
        
        for (int j=0;j<time_step;j++)
        {
            
            double m_sum = accumulate(M_count_statistics[j].begin(), M_count_statistics[j].end(), 0.0);
            double m_mean =  m_sum / init_cell.length(); //按照M_count_statistics的行求均值，即每个cell_collection内所有细胞在特定time_step的均值
            m_means_ratio.pb(m_mean);
            
            double h_sum = accumulate(H_count_statistics[j].begin(), H_count_statistics[j].end(), 0.0);
            double h_mean =  h_sum / init_cell.length();
            h_means_ratio.pb(h_mean);
            
            double u_sum = accumulate(U_count_statistics[j].begin(), U_count_statistics[j].end(), 0.0);
            double u_mean =  u_sum / init_cell.length();
            u_means_ratio.pb(u_mean);
        }
        
        int len_init_cell=init_cell.length();
        int buf_ratio_size=50;
        if ((i >= out_target_start_gen) and (i<= out_target_end_gen))
        {
            for(int t = 0; t < time_step; t++)
            {
                string ratio_file_path=ratio_file_dir+to_string(i)+"_"+to_string(t)+".csv";
                ofstream out_ratio(ratio_file_path,(round_no==round_start)?ios::trunc:ios::app);
                if(!out_ratio){
                    cout << "Unable to open out" << ratio_file_path << endl;
                    exit(1);
                }
                else{
                    char buffer[buf_ratio_size];
                    sprintf(buffer, "%d,%.4f,%.4f,%.4f\n",round_no, m_means_ratio[t],h_means_ratio[t],u_means_ratio[t]);
                    out_ratio<<buffer;
                    out_ratio.close();
                }
            }
            
            int buf_detail_size=len_init_cell+50;
            for(int idx = 0; idx < cell_collection.size();idx++)
            {
                for(int j= 0; j < time_step;j++)
                {
                    string detail_file_path=detail_file_dir+to_string(i)+"_"+to_string(j)+".csv";
                    ofstream out_detail(detail_file_path,(round_no==round_start)?ios::trunc:ios::app);
                    
                    if(!out_detail){
                        cout << "Unable to open out" << detail_file_path << endl;
                        exit(1);
                    }
                    else{
                        char buffer[buf_detail_size];//在index的size不等于1时,得到的结果是多细胞的
                        sprintf(buffer, "%d,%s\n", round_no, out_detail_seq_arr[0][j].c_str());
                        out_detail<<buffer;
                        out_detail.close();
                    }
                }
            }
        }
        
        
        cell_collection = cells_wait_to_add;
        clock_t gen_end_time=clock();
        cout<< "generation time:"<<static_cast<double>(gen_end_time-gen_start_time)/CLOCKS_PER_SEC<<"s"<<endl;//输出1代运行时间
    }
}

vector<int> get_pos_list_from_bed_file(string input_bed_file_path,int max_cpg_sites=1000)
{
    ifstream bed_file(input_bed_file_path);
    int  buf_size=100;
    char buffer[buf_size];
    vector<int> pos_list;
    
    if (!bed_file)
    {
        cout << "Unable to open " << input_bed_file_path << endl;
    }
    else
    {
        int pos=0;
        float methy_level=0.0;
        int cnt=0;
        
        while (! bed_file.eof() )
        {
            if (max_cpg_sites < 0 or (cnt < max_cpg_sites and max_cpg_sites > 0))
            {
                bed_file.getline(buffer,buf_size);
                sscanf(buffer,"%d %f",&pos,&methy_level);
                pos_list.pb(pos);
                cnt++;
            }
            else
            {
                break;
            }
        }
        bed_file.close();
        printf("get pos list from %s successfully!\n",input_bed_file_path.c_str() );
    }
    return pos_list;
}

string generate_CpG_in_methylation_percent_UHM(int CpG_sites_counts,double m_ratio, double u_ratio)
{
    string CpG_str;
    
    for (int i = 0; i < CpG_sites_counts; i++)
    {
        double random_num = (double)(rand()%RAND_MAX)/RAND_MAX;
        
        if (random_num < m_ratio)
        {
            CpG_str.append(1, 'M');
        }
        else if (random_num > m_ratio and random_num < m_ratio + u_ratio)
        {
            CpG_str.append(1, 'U');
        }
        else
        {
            CpG_str.append(1, 'H');
        }
    }
    
    return CpG_str;// 返回生成的状态字符串, 长度=CpG_sites_counts
}

void load_norm_distri_param(string norm_param_file_path)
{
    ifstream param_file(norm_param_file_path);
    int  buf_size=50;
    char buffer[buf_size];
    
    if (!param_file)
    {
        cout << "Unable to open " << norm_param_file_path << endl;
    }
    else
    {
        map<string,float> params_map;
        
        cmatch match_result;  //保存匹配结果
        string pattern="([a-z]+)=([\\d]+\\.[\\d]*)";
        regex regex_expression(pattern);
        while (! param_file.eof() )
        {
            param_file.getline(buffer,buf_size);
            bool match=regex_search(buffer,
                                    buffer+strlen(buffer),
                                    match_result,
                                    regex_expression
                                    );
            if (match)
            {
                string param_name=match_result[1].str();
                float param_value;
                sscanf(match_result[2].str().c_str(),"%f",&param_value);
                params_map[param_name]=param_value;
                //                cout << param_name << " = " << param_value << endl;
            }
        }
        param_file.close();
        
        beta=params_map["beta"];
        gamma_param=params_map["gamma"];
        power=params_map["power"];
        alpha=params_map["alpha"];
        c_off=params_map["c_off"];
        printf("get norm params from %s successfully!\n",norm_param_file_path.c_str() );
    }
    
}

vector<double> load_reaction_param(string reaction_param_file_path)
{
    ifstream param_file(reaction_param_file_path);
    int  buf_size=50;
    char buffer[buf_size];
    vector<double> reaction_params;
    if (!param_file)
    {
        cout << "Unable to open " << reaction_param_file_path << endl;
    }
    else
    {
        map<string,float> params_map;
        
        cmatch match_result;  //保存匹配结果
        string pattern="([A-Za-z_\\+|-]+)=([\\d]+\\.[\\d]*)";
        regex regex_expression(pattern);
        while (! param_file.eof() )
        {
            param_file.getline(buffer,buf_size);
            bool match=regex_search(buffer,
                                    buffer+strlen(buffer),
                                    match_result,
                                    regex_expression
                                    );
            if (match)
            {
                string param_name=match_result[1].str();
                float param_value;
                sscanf(match_result[2].str().c_str(),"%f",&param_value);
                params_map[param_name]=param_value;
                cout << param_name << " = " << param_value << endl;
            }
        }
        param_file.close();
        
        double U_plus_in=params_map["U_plus_in"];
        reaction_params.pb(U_plus_in);
        double H_plus_in=params_map["H_plus_in"];
        reaction_params.pb(H_plus_in);
        double M_minus_in=params_map["M_minus_in"];
        reaction_params.pb(M_minus_in);
        double H_minus_in=params_map["H_minus_in"];
        reaction_params.pb(H_minus_in);
        double H_p_H_in=params_map["H_p_H_in"];
        reaction_params.pb(H_p_H_in);
        double H_p_M_in=params_map["H_p_M_in"];
        reaction_params.pb(H_p_M_in);
        double U_p_M_in=params_map["U_p_M_in"];
        reaction_params.pb(U_p_M_in);
        double H_m_U_in=params_map["H_m_U_in"];
        reaction_params.pb(H_m_U_in);
        double M_m_U_in=params_map["M_m_U_in"];
        reaction_params.pb(M_m_U_in);
        printf("get reaction params from %s successfully!\n",reaction_param_file_path.c_str() );
    }
    return reaction_params;
}

vector<int> construct_n_cpg_sites_for_exp_distribution(int max_cpg_sites, float geometric_p)
{
    vector<int> index_pos_list;
    int first_cpg_pos = rand() % 9 + 1;
    index_pos_list.pb(first_cpg_pos);
    
    //    cout << 0 << " : " << first_cpg_pos <<endl;
    
    default_random_engine generator(time(NULL));
    geometric_distribution<int> distribution(geometric_p);
    int number;
    for (int i=1; i< max_cpg_sites; ++i) {
        number= distribution(generator)+2;
        first_cpg_pos=first_cpg_pos+number;
        //        cout << i << " : " << first_cpg_pos <<endl;
        index_pos_list.pb(first_cpg_pos);
    }
    
    return index_pos_list;
}

void sort_to_bed(int round_size,vector<int> pos_list,string sort_detail_dir,string bed_files_dir,int generation_start,int generation_end,int time_steps,int max_cpg_sites)
{
    int  buf_size=max_cpg_sites+50;
    char buffer[buf_size];
    
    for(int gen=generation_start;gen<=generation_end;gen++)
    {
        for (int t=0;t<time_steps;t++)
        {
            int round=0;
            unsigned long len_methy_seq=0;
            vector<int> methy_status;
            
            string detail_file_path=sort_detail_dir+to_string(gen)+"_"+to_string(t)+".csv";
            ifstream detail_file(detail_file_path);
            
            if(!detail_file){
                cout << "Unable to open read " << detail_file_path << endl;
                exit(1);
            }
            else{
                string bed_file_path=bed_files_dir+to_string(gen)+"_"+to_string(t)+".bed";
                ofstream bed_file(bed_file_path);
                
                char methy_seq[buf_size];
                int line_cnt=0;
                while (! detail_file.eof() && line_cnt< round_size)
                {
                    line_cnt=line_cnt+1;
                    detail_file.getline(buffer,buf_size);
                    sscanf(buffer,"%d,%s\n",&round,methy_seq);
                    len_methy_seq=strlen(methy_seq);
                    for(unsigned long i=0;i<len_methy_seq;i++)
                    {
                        char ch=methy_seq[i];
                        if (ch=='M')
                        {
                            if(line_cnt==1)
                            {
                                methy_status.push_back(2);
                            }
                            else
                            {
                                methy_status[i]=methy_status[i]+2;
                            }
                        }
                        else if (ch=='H')
                        {
                            if(line_cnt==1)
                            {
                                methy_status.push_back(1);
                            }
                            else
                            {
                                methy_status[i]=methy_status[i]+1;
                            }
                        }
                        else if (ch=='U')
                        {
                            if(line_cnt==1)
                            {
                                methy_status.push_back(0);
                            }
                        }
                    }
                }
//                vector<float> methy_mean;
                char wrt_buffer[100];
                float round_cnt=2*(float)round;
                for(unsigned long i=0;i<len_methy_seq;i++)
                {
//                    methy_mean.pb(methy_status[i]/round_cnt);
                    sprintf(wrt_buffer,"%d %.6f\n",pos_list[i],methy_status[i]/round_cnt);
                    bed_file<<wrt_buffer;
                }
                detail_file.close();
                bed_file.close();
            }
            printf( "finished convert %d.%d to bed\n",gen,t);
        }
        
    }
}

map<int,float> read_bed_file_and_store_pos_to_a_struct(string bed_file_path,bool ignore_d)
{
    map<int,float> struct_to_store;
    char buffer[100];
    int pos;
    float methy_level;
    int index=0;
    ifstream bed_file(bed_file_path);
    if(!bed_file){
        cout << "Unable to open read " << bed_file_path << endl;
        exit(1);
    }
    else{
        while (!bed_file.eof() )
        {
            index=index+1;
            bed_file.getline(buffer,100);
            sscanf(buffer,"%d %f",&pos,&methy_level);
            if (ignore_d)
            {
                pos=index;
            }
            struct_to_store[pos]=methy_level;
        }
        bed_file.close();
    }
    return struct_to_store;
}

vector<vector<float>> filter_d_length_to_generate_CpG_pairs_not_inter_with_other_cpg(int d,vector<int> keys,vector<float> vals)
{
    int key_size=(int)keys.size();
    int pre_key,post_key;
    float pre_val,post_val;
    vector<vector<float>> array_to_store_pairs(2);
    
    for (int i=0;i< key_size;i++)
    {
        pre_key = keys[i];
        post_key = keys[i+1];
        
        if (pre_key + d == post_key)
        {
            pre_val=vals[i];
            post_val=vals[i+1];
            //            printf("<%d,%d> : %.2f, %.2f\n",pre_key,post_key,pre_val,post_val);
            array_to_store_pairs[0].push_back(pre_val);
            array_to_store_pairs[1].push_back(post_val);
        }
    }
    
    return array_to_store_pairs;
}
vector<vector<float>> filter_d_length_to_generate_CpG_pairs(int d,map<int,float> CpG_pos_and_methy_struct)
{
    int map_size=(int)CpG_pos_and_methy_struct.size();
    int pre_key,post_key;
    float pre_val,post_val;
    vector<vector<float>> array_to_store_pairs(2);
    for(map<int,float>::iterator it = CpG_pos_and_methy_struct.begin(); it != CpG_pos_and_methy_struct.end(); ++it) {
        pre_key=it->first;
        pre_val=it->second;
        post_key=pre_key+d;
        if (CpG_pos_and_methy_struct.find(post_key) != CpG_pos_and_methy_struct.end())
        {
            post_val=CpG_pos_and_methy_struct[post_key];
            array_to_store_pairs[0].push_back(pre_val);
            array_to_store_pairs[1].push_back(post_val);
        }
    }
    return array_to_store_pairs;
}

float calc_C_d_by_pearson_correlation(vector<vector<float>> CpG_pairs)
{
    float sum_pre=0.0;
    float sum_post=0.0;
    
    int length = (int)CpG_pairs[0].size();
    for (int i=0;i<length;i++)
    {
        sum_pre=sum_pre+CpG_pairs[0][i];
        sum_post=sum_post+CpG_pairs[1][i];
    }
    
    float mean1=sum_pre/float(length);
    float mean2=sum_post/float(length);
    
    float sum_up=0.0;
    float sum_down_left=0.0;
    float sum_down_right=0.0;
    
    float xi,yi;
    for (int i=0;i<length;i++)
    {
        xi=CpG_pairs[0][i];
        yi=CpG_pairs[1][i];
        
        sum_up=sum_up+(xi-mean1)*(yi-mean2);
        sum_down_left=sum_down_left+(xi-mean1)*(xi-mean1);
        sum_down_right=sum_down_right+(yi-mean2)*(yi-mean2);
    }
    
    float sum_down=sqrt(sum_down_left*sum_down_right);
    
    if((sum_down >= - EPSINON) && (sum_down <= EPSINON))
    {
        return -5;
    }
    float rd=sum_up/sum_down;
    return rd;
}

void calc_correlation(string bed_file_path,string rd_file_path,int d_max,bool is_inter_with_other_cpg,bool ignore_d=false)
{
    ofstream rd_file(rd_file_path);
    
    if (!rd_file)
    {
        cout << "Unable to open read " << rd_file_path << endl;
        exit(1);
    }
    
    map<int,float> CpG_pos_and_methy_struct = read_bed_file_and_store_pos_to_a_struct(bed_file_path, ignore_d);
    
    vector<int> keys;
    vector<float> vals;
    for(map<int,float>::iterator it = CpG_pos_and_methy_struct.begin(); it != CpG_pos_and_methy_struct.end(); ++it) {
        keys.push_back(it->first);
        vals.push_back(it->second);
        //        cout << it->first<<" "<< it->second << endl;
    }
    
    vector<vector<float>> CpG_pairs;
    int d_count=0;
    float rd=0.0;
    char ltw[100];
    for (int d=2;d<d_max;d++)
    {
        if (is_inter_with_other_cpg)
        {
            CpG_pairs=filter_d_length_to_generate_CpG_pairs(d,CpG_pos_and_methy_struct);
        }
        else
        {
            CpG_pairs=filter_d_length_to_generate_CpG_pairs_not_inter_with_other_cpg(d,keys,vals);
        }
        d_count = (int)CpG_pairs.size();
        
        if (d_count)
        {
            rd=calc_C_d_by_pearson_correlation(CpG_pairs);
            if (rd > -2)
            {
                printf("chr%d d=%d, rd=%f\n",1,d,rd);
                sprintf(ltw,"%d,%f\n",d,rd);
                rd_file<<ltw;
            }
        }
        else{
            printf("chr%d passed d=%d",1,d);
        }
    }
    rd_file.close();
    
    printf("\n\n");
}

void calc_mean_rd_from_rd_dir(string rd_dir_name,string out_file_path,int d_max,int generation_start,int generation_end,int time_steps)
{
    
    string rd_file_path;
    ofstream out_file(out_file_path);
    if (!out_file)
    {
        cout << "Unable to open out " << out_file_path << endl;
        exit(1);
    }
    vector<vector<float>> rd_vect;
    for (int i=0;i<d_max;i++)
    {
        vector<float> tmp_vect;
        rd_vect.push_back(tmp_vect);
    }
    
    int d=-1;
    float methy_level=0.0;
    char buffer[50];
    float m_sum,m_mean;
    for(int gen=generation_start;gen <= generation_end;gen++)
    {
        for(int step = 0; step < time_steps; step++)
        {
            rd_file_path=rd_dir_name+to_string(gen)+"_"+to_string(step)+".csv";
            ifstream rd_file(rd_file_path);
            if (!rd_file)
            {
                cout << "Unable to open read " << rd_file_path << endl;
                exit(1);
            }
            while (!rd_file.eof() )
            {
                rd_file.getline(buffer,50);
                sscanf(buffer,"%d,%f",&d,&methy_level);
                rd_vect[d].pb(methy_level);
            }
            rd_file.close();
        }
    }
    vector<float> rd_list;
    char ltw[50];
    for (d=2;d<d_max;d++)
    {
        m_sum = accumulate(rd_vect[d].begin(), rd_vect[d].end(), 0.0);
        m_mean =  m_sum / (float) rd_vect[d].size(); //按照M_count_statistics的行求均值，即每个cell_collection内所有细胞在特定time_step的均值
        sprintf(ltw,"%d,%f\n",d,m_mean);
        out_file<<ltw;
    }
    out_file.close();
    printf("writing mean rd finished!\n");
}

void calc_correlation_for_generations(int generation_start,int generation_end,int time_steps,string bed_dir,string rd_without_dir,int d_max,bool calc_interval=false,bool ignore_d=false)
{
    for (int gen=generation_start; gen <= generation_end; gen++)
    {
        for(int t=0; t< time_steps;t++)
        {
            string input_bed_file_path=bed_dir+to_string(gen)+"_"+to_string(t)+".bed";
            string out_rd_file_path=rd_without_dir+to_string(gen)+"_"+to_string(t)+".csv";
            calc_correlation(input_bed_file_path, out_rd_file_path, d_max, calc_interval,ignore_d);
        }
    }
    
}
void start_simulation()
{
    int repeat_start=1;
    int repeat_end=1;
    
    int generations=30;
    int round_start=1;
    int round_end=6;
    int round_size=round_end-round_start+1;
    int max_cpg_sites=100000;
    string init_cell;
    double m_ratio=0.181214;
    double u_ratio=0.391004;
    vector<int> index_pos_list;
    string path_dir="/Users/Ren/XCodeProjects/Methylation_Server/Methylation_Server/";
    string input_bed_file_path=path_dir+"chr1.bed";
    string ratio_file_dir;
    string detail_file_dir;
    string bed_file_dir;
    string rd_without_dir;
    int time_step=100;
    int d_max=1000;
    bool calc_interval=false;
    bool ignore_d=false;
    load_norm_distri_param("/Users/Ren/PycharmProjects/Methylation_Server/input_new/diff_period_param/pat_2.txt");
    vector<double> propensity_list=load_reaction_param("/Users/Ren/PycharmProjects/Methylation_Server/input_new/reaction_0.txt");
    int nearby = 1;
    int max_cells = 2;
    bool real_nearby = false;
    
    bool simulation=true;
    bool calc_corr=true;
    string out_mean_rd_file;
    int out_target_start_gen=generations;
    int out_target_end_gen=generations;
    
    for (int rp=repeat_start;rp<=repeat_end;rp++){
        int partial_start=1;
        int partial_end=1;
        path_dir=path_dir+"repeat_"+to_string(rp)+"/";
        int rtn_code=exist_or_make_dir(path_dir);
        for (int partial_i=partial_start;partial_i<=partial_end;partial_i++){
            path_dir=path_dir+"partial_"+to_string(partial_i)+"/";
            rtn_code=exist_or_make_dir(path_dir);
            bool real_chr_pos=true;
            if(!real_chr_pos){
                //超几何分布参数
                float geometric_p = 0.3;
                index_pos_list=construct_n_cpg_sites_for_exp_distribution(max_cpg_sites, geometric_p);
            }
            else{
                index_pos_list=get_pos_list_from_bed_file(input_bed_file_path,max_cpg_sites);
            }
            
            init_cell=generate_CpG_in_methylation_percent_UHM(max_cpg_sites,m_ratio,u_ratio);
            ratio_file_dir=path_dir+"ratio/";//detail的生成序列
            rtn_code=exist_or_make_dir(ratio_file_dir);
            detail_file_dir=path_dir+"detail/";//detail的生成序列
            rtn_code=exist_or_make_dir(detail_file_dir);
            
            if (simulation){
                for(int round_i=round_start;round_i<=round_end;round_i++)
                {
                    simulate(round_i,round_start,out_target_start_gen,out_target_end_gen,generations,time_step,init_cell,detail_file_dir,ratio_file_dir,propensity_list,index_pos_list,nearby,max_cells,real_nearby);
                }
            }
            if(calc_corr){
                bed_file_dir=path_dir+"bed/";//detail的生成序列
                rtn_code=exist_or_make_dir(bed_file_dir);
                
                rd_without_dir=path_dir+"rd_without/";
                rtn_code=exist_or_make_dir(rd_without_dir);
                
                sort_to_bed(round_size,index_pos_list,detail_file_dir,bed_file_dir,out_target_start_gen,out_target_end_gen,time_step,max_cpg_sites);
                calc_correlation_for_generations(out_target_start_gen,out_target_end_gen,time_step,bed_file_dir,rd_without_dir,d_max,calc_interval,ignore_d=false);
                out_mean_rd_file=path_dir+"rd_mean.csv";
                calc_mean_rd_from_rd_dir(rd_without_dir,out_mean_rd_file,d_max,out_target_start_gen,out_target_end_gen,time_step);
            }
            
        }
    }
}


int main(int argc, const char * argv[]) {
    start_simulation();
    return 0;
}
