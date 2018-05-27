
# coding: utf-8

# In[3]:


from Bio import SeqIO 
import pandas as pd

class kmer:
    counter = 1 
    sequence = ''
    locuses = []
    end_locuses =[]
    
    def __init__(self, kmer_name): # последовательность кмера
        self.sequence = kmer_name
        
    def  increase(self): #счетчик
        self.counter += 1
        
    def increase_n(self, n): #инкрементация счетчика
        self.counter += n
    
    def find_locuses(self, seq): #поиск локусов (начальных и конечных координат)
        l = len(self.sequence)
        for index in range(len(seq)-l+1):
            current_kmer = seq[index:(index+l)]
            if current_kmer == self.sequence:
                self.locuses.append(index)
                self.end_locuses.append(index+l)
    
    def output(self): #вывод локусов в виде дата фрэйма
        out_table = pd.DataFrame(data = {'start_positions': self.locuses, 'end_positions': self.end_locuses})
        return (out_table)

    
seq = ''

with open ("/Users/yukornienko/Downloads/seq_y_pestis.fasta") as fasta:
    for record in SeqIO.parse(fasta, "fasta"):
        seq = record.seq
        seq_len = len(seq)
        kmer_size = 23
        kmer_dict = {}
#сохраняем текущую хромосому
        chrom = record.name
        
        for index in range(seq_len-kmer_size+1): #создаем словарь кмеров с числом их встречания в последовательности
            current_kmer = seq[index:(index+kmer_size)]
            if current_kmer in kmer_dict:
                kmer_dict[current_kmer].increase()
            else:
                kmer_dict[current_kmer] = kmer(current_kmer)
                
        most_freq_kmer_seq= max(kmer_dict.keys(), key=(lambda key: kmer_dict[key].counter)) #самый частый кмер
        #самый частый кмер - вывод его последовательности, числа раз встречания и координат
        
        #вывод вместе с хромосомой
        print(chrom, most_freq_kmer_seq)
        print(kmer_dict[most_freq_kmer_seq].counter)
        
        most_freq_kmer = kmer(most_freq_kmer_seq)
        most_freq_kmer.find_locuses(seq)
       
        
        most_freq_kmer.output()

