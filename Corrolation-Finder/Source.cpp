/*
* Core.cpp
*
*  Created on: 1 Feb. 2018
*      Author: Arash Bayat
*/

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define Error(X) {printf("*** ERROR: %s (line:%u - File %s)\n", X, __LINE__, __FILE__); exit(0);}
#define NULL_CHECK(X) {if(!X) {printf("*** ERROR: %s is null (line:%u - File %s)\n", #X, __LINE__, __FILE__); exit(0);}}
#define Printf printf
//#define PrintfD printf
#define PrintfD DO_NOTHING

void DO_NOTHING(const char *x, ...)
{
	return;
}

typedef char SAMPLE_CLASS;		// the data type to hold sample class label
typedef unsigned int	NODE_ID;		// the data type to keep node id.
typedef char VAR_VALUE;			// the data type to keep the value of each variable
typedef char VAR_ORD;			// the data type to keep ordinality of variables.
typedef unsigned int SAMPLE_COUNTER;		// the data type to keep number of samples
typedef unsigned int uint;
typedef unsigned long long int ulint;

typedef enum
{
	memory,					// variable data stored in a table in the memory
	compressed_memory,		// variable data table is compressed in the memory
	RASM						// random access to records with sequential access to file
} MEMORY_PLAN;

#define MAX_FILE_NAME_LEN 4096

#define MAX_NAME_LEN 32
typedef struct
{
	char name[MAX_NAME_LEN];
} NAME;

typedef struct
{
	double InfoGaind;		// InfoGaind independently by this value
	double sumInfoGained;	// total amount of InfoGained by this variable in tree
	uint numInfoGained;		// number of time this varibale is used in forest.
	void print(char *varName)
	{
		printf("\n%s\t%20.19f\t%20.19f\t%u\t%20.19f", varName, InfoGaind, sumInfoGained/numInfoGained, numInfoGained, sumInfoGained);
	}

} VAR_DATA;

ulint FILE_SIZE(char *fileName)
{
	FILE *f = fopen(fileName, "r");
	NULL_CHECK(f);
	fseek(f, 0L, SEEK_END);
	ulint size = ftell(f);
	fclose(f);
	return size;
}

ulint FILE_LINE(char *fileName)
{
	FILE *f = fopen(fileName, "r");
	NULL_CHECK(f);
	int ch = 0;
	ulint lines = 0;

	lines++;
	while ((ch = fgetc(f)) != EOF)
	{
		if (ch == '\n')
			lines++;
	}
	fclose(f);
	return lines;
}

typedef struct T_UNCOMPRESSED_MEMORY
{
	SAMPLE_COUNTER num_sample;	// number of samples
	uint num_var;				// number of variables

	VAR_VALUE *data;		// keep value for each variable (it is a 2D array see Get_Var_Value())
	VAR_ORD *ord;			// ordinallity of each variable
	SAMPLE_CLASS *sample_class; // class labels for samples

	NAME *varName;
	NAME *sampleName;

	uint ord_class;				// Ordinality of classes. How many class exist in total.
	double sample_purity;		// purity of samples
	
	void ComputeSamplePurity()
	{
		SAMPLE_COUNTER *cnt = new SAMPLE_COUNTER[ord_class];
		memset(cnt, 0, sizeof(SAMPLE_COUNTER) * ord_class);
		for (uint i = 0; i < num_sample; i++)
		{
			cnt[sample_class[i]]++;
		}

		sample_purity = 0;
		for (uint i = 0; i<ord_class; i++)
		{
			double fij = (double)cnt[i] / num_sample; // fraction of sample in node i with class j;
			sample_purity += (fij * fij);
		}
	}

	VAR_VALUE *GetVarValue(uint index)
	{
		return &data[index * num_sample];
	}

	void FillRandom(uint a_num_var, SAMPLE_COUNTER a_num_samples, VAR_ORD a_ord)
	{
		num_var = a_num_var;
		num_sample = a_num_samples;

		// allocate memory
		data = new VAR_VALUE[num_var * num_sample];
		NULL_CHECK(data);
		ord = new VAR_ORD[num_var];
		NULL_CHECK(ord);
		sample_class = new SAMPLE_CLASS[num_sample];
		NULL_CHECK(sample_class);
		sampleName = new NAME[num_sample];
		NULL_CHECK(sampleName);
		varName = new NAME[num_var];
		NULL_CHECK(varName);

		// fill random value for variables
		for (uint i = 0; i<num_var; i++)
		{
			ord[i] = a_ord;
			for (uint j = 0; j<num_sample; j++)
				data[(i * num_sample) + j] = rand() % ord[i];
			sprintf(varName[i].name, "VAR_%010u", i);
		}

		ord_class = 2;
		// assign random class for samples
		for (uint i = 0; i<num_sample; i++)
		{
			sample_class[i] = rand() % 2;
			sprintf(sampleName[i].name, "SAMPEL_%010u", i);
		}
	}

	// this function fill read VCF data variables are one byte
	void ParseTable(char *prefix)
	{
		if (sizeof(VAR_VALUE) != 1)
		{
			Error("This function only works when var value type is one byte");
		}

		char fileName[MAX_FILE_NAME_LEN];

		// count number of samples and read sample names
		sprintf(fileName, "%s.Sample", prefix);
		num_sample = (SAMPLE_COUNTER) FILE_LINE(fileName) - 1;
		sampleName = new NAME[num_sample];
		NULL_CHECK(sampleName);
		Printf("\nThere are %u samples in %s", num_sample, fileName);
		FILE *fsam = fopen(fileName, "r");
		NULL_CHECK(fsam);
		for (uint i = 0; i < num_sample; i++)
		{
			fscanf(fsam, "%s", sampleName[i].name);
		}
		fclose(fsam);

		// count number of sites and read site names
		sprintf(fileName, "%s.Site", prefix);
		num_var = (uint) FILE_LINE(fileName) - 1;
		varName = new NAME[num_var];
		NULL_CHECK(varName);
		Printf("\nThere are %u sites in %s", num_var, fileName);
		FILE *fvar = fopen(fileName, "r");
		NULL_CHECK(fvar);
		for (uint i = 0; i < num_var; i++)
		{
			fscanf(fvar, "%s", varName[i].name);
		}
		fclose(fvar);

		ulint gt = (ulint)num_sample * num_var; // number of genotype should exist in Genotype file
		
		sprintf(fileName, "%s.GT", prefix);
		ulint num_gt = FILE_SIZE(fileName);
		Printf("\nThere are %llu genotypes in %s (%u*%u=%llu)", num_gt, fileName, num_sample, num_var, gt);
		if (gt != num_gt)
		{
			Error("Number of genotypes does not match number of samples and variants");
		}
		
		// read genotype file
		data = new VAR_VALUE[num_gt];
		NULL_CHECK(data);
		FILE *fgt = fopen(fileName, "rb");
		NULL_CHECK(fgt);
		fread(data, 1, num_gt, fgt);
		// can be improved by SSE
		for (ulint i = 0; i < num_gt; i++)
		{
			data[i] -= '0';
		}

		ord = new VAR_ORD[num_var];
		NULL_CHECK(ord);
		for (uint i = 0; i < num_var; i++)
		{
			ord[i] = 3;
		}

		return;
	}

	void SimLabelNoCorrolation(uint num_rand_site, char *prefix)
	{
		ord_class = 2;

		sample_class = new SAMPLE_CLASS[num_sample];
		NULL_CHECK(sample_class);

		uint *rnd = new uint[num_rand_site];
		NULL_CHECK(rnd);

		// select variant for simulation
		for (uint i = 0; i < num_rand_site; i++)
			rnd[i] = rand() % num_var;

		uint num_labeled = 0;

		for (uint i = 0; i < num_sample; i++)
		{
			sample_class[i] = 0;
			for(uint j=0; j<num_rand_site; j++)
				if (data[(rnd[j] * num_sample) + i] != 0) // if it is not homoref
				{
					sample_class[i] = 1;
					num_labeled++;
					break;
				}
		}
		
		Printf("\nNumber of labled sample: %u", num_labeled);

		char fileName[MAX_FILE_NAME_LEN];

		// writing labels
		sprintf(fileName, "%s.sim.csv", prefix);
		FILE *fsim = fopen(fileName, "w");
		fprintf(fsim, "Sample,label\n");
		for (uint i = 0; i < num_sample; i++)
		{
			fprintf(fsim, "%s,%u\n", sampleName[i].name, sample_class[i]);
		}
		fclose(fsim);

		// writing selected variants
		sprintf(fileName, "%s.sim.var", prefix);
		FILE *fvar = fopen(fileName, "w");\
		fprintf(fvar, "index\tCHR_POS_ALT\n");
		for (uint i = 0; i < num_rand_site; i++)
		{
			fprintf(fvar, "%u\t%s\n", rnd[i], varName[rnd[i]].name);
		}
		fclose(fvar);
		
		delete[] rnd;
	}

	void SimLabelCorrolation(uint num_rand_site, char *prefix)
	{
		ord_class = 2;

		sample_class = new SAMPLE_CLASS[num_sample];
		NULL_CHECK(sample_class);

		uint *rnd = new uint[num_rand_site];
		NULL_CHECK(rnd);

		// select variant for simulation
		for (uint i = 0; i < num_rand_site; i++)
			rnd[i] = rand() % num_var;

		uint num_labeled = 0;

		for (uint i = 0; i < num_sample; i++)
		{
			sample_class[i] = 1; // assume the lable is one. then we set it as 0 if any of the condition does not work
			for (uint j = 0; j<num_rand_site; j++)
			{
				VAR_VALUE gt = data[(rnd[j] * num_sample) + i]; // the genotype
				if (j & 1) // odd number
				{
					if (gt != 0)
						sample_class[i] = 0;
				}
				else
				{
					if (gt == 0)
						sample_class[i] = 0;
				}
			}
			if (sample_class[i] == 1)
				num_labeled++;
		}

		Printf("\nNumber of labled sample: %u", num_labeled);

		char fileName[MAX_FILE_NAME_LEN];

		// writing labels
		sprintf(fileName, "%s.sim.csv", prefix);
		FILE *fsim = fopen(fileName, "w");
		fprintf(fsim, "Sample,label\n");
		for (uint i = 0; i < num_sample; i++)
		{
			fprintf(fsim, "%s,%u\n", sampleName[i].name, sample_class[i]);
		}
		fclose(fsim);

		// writing selected variants
		sprintf(fileName, "%s.sim.var", prefix);
		FILE *fvar = fopen(fileName, "w"); \
			fprintf(fvar, "index\tCHR_POS_ALT\n");
		for (uint i = 0; i < num_rand_site; i++)
		{
			fprintf(fvar, "%u\t%s\n", rnd[i], varName[rnd[i]].name);
		}
		fclose(fvar);

		delete[] rnd;
	}

	void Loadlabel(char *fileName)
	{
		ord_class = 2;

		sample_class = new SAMPLE_CLASS[num_sample];
		NULL_CHECK(sample_class);

		// read labels
		FILE *flbl = fopen(fileName, "r");

		for (uint i = 0; i < num_sample; i++)
		{
			int temp;
			fscanf(flbl, "%u", &temp);
			sample_class[i] = temp;
		}
		fclose(flbl);

		return;
	}

	void StoreToFile(const char *prefix)
	{
		char fn[MAX_FILE_NAME_LEN];
		sprintf(fn, "%s.vds", prefix);
		FILE *file = fopen(fn, "wb");
		NULL_CHECK(file);
		fwrite(this, sizeof(T_UNCOMPRESSED_MEMORY), 1, file);
		fwrite(this->data, sizeof(VAR_VALUE), (num_var * num_sample), file);
		fwrite(this->ord, sizeof(VAR_ORD), (num_var), file);
		fwrite(this->sample_class, sizeof(SAMPLE_CLASS), (num_sample), file);
		fwrite(this->sampleName, sizeof(NAME), (num_sample), file);
		fwrite(this->varName, sizeof(NAME), (num_var), file);
	}

	void LoadFromFile(const char *prefix)
	{
		char fn[MAX_FILE_NAME_LEN];
		sprintf(fn, "%s.vds", prefix);
		FILE *file = fopen(fn, "rb");
		NULL_CHECK(file);
		fread(this, sizeof(T_UNCOMPRESSED_MEMORY), 1, file);

		data = new VAR_VALUE[num_var * num_sample];
		NULL_CHECK(data);
		ord = new VAR_ORD[num_var];
		NULL_CHECK(ord);
		sample_class = new SAMPLE_CLASS[num_sample];
		NULL_CHECK(sample_class);
		sampleName = new NAME[num_sample];
		NULL_CHECK(sampleName);
		varName = new NAME[num_var];
		NULL_CHECK(varName);

		fread(this->data, sizeof(VAR_VALUE), (num_var * num_sample), file);
		fread(this->ord, sizeof(VAR_ORD), (num_var), file);
		fread(this->sample_class, sizeof(SAMPLE_CLASS), (num_sample), file);
		fread(this->sampleName, sizeof(NAME), (num_sample), file);
		fread(this->varName, sizeof(NAME), (num_var), file);
	}

	~T_UNCOMPRESSED_MEMORY()
	{
		delete[] data;
		delete[] ord;
		delete[] sample_class;
		delete[] varName;
		delete[] sampleName;
	}
} UNCOMPRESSED_MEMORY;

typedef struct
{
	uint ord_var;	// tree[i].ord_var=j indicates that variable used in depth 'i' to split data has the ordinality of 'j'. Thus the value of the variable is between 0 and (ord_var[i]-1)
	uint var_id;	// tree[i].var_id=j indicates that id of variable used in depth 'i' to split samples is 'j'.
	double purity;	// tree[i].purity=j indicates that the purity of tree in depth 'i' is 'j'.
} TREE;

class CORE
{
public:

	CORE()
	{
		num_sample = 0;
		sample_class = NULL;
		ord_class = 0;
		max_num_node = 0;
		num_node = 0;
		node_id = NULL;
		tbl_idx = NULL;
		var_value = NULL;
		table = NULL;
		max_depth = 0;
		depth = 0;
		tree = NULL;
		num_tree_per_file = 0;
		num_stored_tree = 0;
		file_tree = NULL;
		tree_file_name_prefix[0] = 0;
		file_counter = 0;
		thread_id = -1;
		num_var = 0;
		memory_plan = memory;
		uncompressed_memory = NULL;
		varData = NULL;
		sample_purity = 0;
		varName = NULL;
		sampleName = NULL;
	}

	~CORE()
	{
		if (!node_id)
			delete[] node_id;
		if (!tbl_idx)
			delete[] tbl_idx;
		if (!table)
			delete[] table;
		if (!tree)
			delete[] tree;
		if (!varData)
			delete[] varData;

		if (file_tree) // not to leave output file open
			fclose(file_tree);
	}

	void Init(SAMPLE_COUNTER a_num_sample, SAMPLE_CLASS *a_sample_class, uint a_ord_class, double a_sample_purity, uint a_max_depth, uint a_num_tree_per_file, uint a_num_var, uint a_max_var_ord, MEMORY_PLAN a_memory_plan, void *a_dataset, char *prefix, uint a_thread_id)
	{
		num_sample = a_num_sample;
		sample_class = a_sample_class;

		ord_class = a_ord_class;
		max_depth = a_max_depth;

		// TO BE IMPROVED: if there is only one variable with high ordinality the chance of choosing that variable is low. But we assume it is selected in all tree nodes.
		max_num_node = 1;
		for (uint i = 0; i<max_depth; i++)
			max_num_node *= a_max_var_ord;

		node_id = new NODE_ID[num_sample];
		NULL_CHECK(node_id);

		tbl_idx = new NODE_ID[num_sample];
		NULL_CHECK(node_id);

		table = new SAMPLE_COUNTER[max_num_node * ord_class];
		NULL_CHECK(table);

		tree = new TREE[max_depth];
		NULL_CHECK(tree);

		num_tree_per_file = a_num_tree_per_file;
		num_stored_tree = num_tree_per_file; // this is to force opening the first file when storing a tree in file.

		if (strlen(prefix) >= MAX_FILE_NAME_LEN)
			Error("Tree prefix file name is too long. Use shorter path (max 4095 characters)");
		strcpy(tree_file_name_prefix, prefix);
		thread_id = a_thread_id;

		num_var = a_num_var;

		varData = new VAR_DATA[num_var];
		NULL_CHECK(varData);
		memset(varData, 0, sizeof(VAR_DATA) * num_var);

		sample_purity = a_sample_purity;

		memory_plan = a_memory_plan;
		switch (memory_plan)
		{
		case memory:
			uncompressed_memory = (UNCOMPRESSED_MEMORY *)a_dataset;
			varName = uncompressed_memory->varName;
			sampleName = uncompressed_memory->sampleName;
			break;
		case compressed_memory:
			break;
		case RASM:
			break;
		}
		return;
	}

	void BuildForest(uint num_tree)
	{
		Printf("\nCompute purity of all variable ");
		ComputeAllInfoGained();
		
		Printf("\nBuilding trees ");
		for (uint i = 0; i < num_tree; i++)
		{
			BuildTree();
			StoreTree();
		}

		for (uint i = 0; i<num_var; i++)
		{
			//sprintf(varData[i], "VAR_%08u", i);
			if(varData[i].numInfoGained)
				varData[i].print(varName[i].name);
		}
	}

	void ComputeAndPrintAllInfoGained()
	{
		// Compute purity
		ComputeAllInfoGained();
		// Print purity
		for (uint i = 0; i<num_var; i++)
		{
			printf("\n%u\t%10.9f", i, varData[i].InfoGaind);
		}
		return;
	}

private:
	// variables to keep sample info
	SAMPLE_COUNTER num_sample;	// total number of sample in data set.
	SAMPLE_CLASS *sample_class;	// sample_class[i]=j indicates that sample 'i' has the class label 'j'
	uint ord_class;				// Ordinality of classes. How many class exist in total.
								// variables to build the tree
	NODE_ID max_num_node;	// maximum number of node that could exist on one level of tree (nodes at the same depth in the tree). The memory for the table is allocated base on this value. If max depth of tree is N and max ordinality of a variable is M then the suggested value is M^N.
	NODE_ID num_node;		// number of nodes at current depth of tree.
	NODE_ID *node_id;		// node_id[i]=j indicates that sample 'i' is in node 'j'
	NODE_ID *tbl_idx;		// since we keep table in an array we compute table index before accessing the table. This result in better optimisation with SSE instructions
	VAR_VALUE *var_value;	// var_value[i]=j indicates that for sample 'i' the value of the selected variable is 'j'
	SAMPLE_COUNTER *table;	// table[(i*ord_calss)+j]=k indicates that there are 'k' samples labeled with class 'j' in the node with id 'i'. The table that holds the number of each class of samples in each node.

	// variables to keep track of tree
	uint max_depth;	// maximum depth of tree. The memory for the arrays to hold variable id and ordinality at each depth is allocated based on this value.
	uint depth;		// the current depth of tree.
	TREE *tree;

	// variables related to how to store tree in files
	uint num_tree_per_file;  // number of tree to be store in each file.
	uint num_stored_tree;	// number of tree currently stored in the current file.
	FILE *file_tree;			// the file to store tree to
	char tree_file_name_prefix[4096];
	uint file_counter;
	uint thread_id;

	// variables to read values form input data set
	uint num_var;				// total number of variable in data set
	MEMORY_PLAN memory_plan;
	UNCOMPRESSED_MEMORY *uncompressed_memory;

	NAME *varName;
	NAME *sampleName;

	double sample_purity;	// the purity of samples originally

	VAR_DATA *varData;

	void TakeVariable(uint idx)
	{
		// based on memory plan load variable data
		switch (memory_plan)
		{
		case memory:
			tree[depth].var_id = idx;
			tree[depth].ord_var = uncompressed_memory->ord[idx];
			var_value = uncompressed_memory->GetVarValue(idx);
			break;
		case compressed_memory:
			break;
		case RASM:
			break;
		}

		return;
	}

	void TakeRandomVariable()
	{
		// select a random variable
		uint rnd = (uint)(rand());
		rnd *= (uint)(rand());
		rnd %= num_var;

		// take the random varaible
		TakeVariable(rnd);

		return;
	}

	void StoreTree()
	{
		// *** the initial value of num_stored_tree is set to num_tree_per_file so that this function open the first file when the program started.
		// if the file does not reach its capacity
		if (num_stored_tree >= num_tree_per_file)
		{
			// close the current file and ask for a new file. Also clear the counter
			if (file_tree)	// this check is for the first file to open.
				fclose(file_tree);

			// reset counter
			num_stored_tree = 0;

			//opening new file
			char file_name[5000];
			sprintf(file_name, "%s_%04u_%09u.tree.bin", tree_file_name_prefix, thread_id, file_counter);
			file_tree = fopen(file_name, "wb");
			NULL_CHECK(file_tree);

			file_counter++;
		}

		// store the tree in file
		fwrite(tree, sizeof(TREE), max_depth, file_tree);
		num_stored_tree++;

		return;
	}

	void ComputeAllInfoGained()
	{
		// for every variable
		for (uint i = 0; i<num_var; i++)
		{
			memset(tree, 0, sizeof(TREE) * 2);
			memset(node_id, 0, sizeof(NODE_ID) * num_sample);
			depth = 0;
			num_node = 1;
			// take a variable and split samples
			TakeVariable(i);
			SplitAndComputePurity();
			varData[i].InfoGaind = tree[0].purity - sample_purity;
		}
		return;
	}

	void BuildTree()
	{
		memset(tree, 0, sizeof(TREE) * max_depth);
		memset(node_id, 0, sizeof(NODE_ID) * num_sample);
		depth = 0;
		num_node = 1;
		// for each depth of tree
		for (uint i = 0; i<max_depth; i++)
		{
			// take a variable and split samples
			TakeRandomVariable();
			SplitAndComputePurity();
			// keep track of all information gained by a variable
			double InfoGained;
			if (depth == 0)
			{
				InfoGained = tree[depth].purity - sample_purity;
			}
			else
			{
				InfoGained = tree[depth].purity - tree[depth - 1].purity;
			}

			varData[tree[depth].var_id].sumInfoGained += InfoGained; // (-1 * power(InfoGained, 5));
			varData[tree[depth].var_id].numInfoGained++;

			// increase the depth of the tree
			depth++;
		}
		// print the purity of tree
		//printf("\n");
		//for (uint i = 0; i < max_depth; i++)
		//	printf("%5.3f\t", tree[i].purity);
		return;
	}

	// split samples using the current selected variable.
	void SplitAndComputePurity()
	{
		NODE_ID num_node_next = num_node * tree[depth].ord_var; // number of nodes after splitting samples by the selected variable

		// clear table to count samples split by the selected variable
		memset((char *)table, 0, (num_node_next * ord_class * sizeof(SAMPLE_COUNTER)));

		// TO BE IMPROVED using SSE instruction
		for (uint i = 0; i<num_sample; i++)
		{
			// for all samples compute their new node id based on their current node id and the value of the selected variable
			node_id[i] = ((node_id[i] * tree[depth].ord_var) + var_value[i]);
		}

		// TO BE IMPROVED using SSE instruction
		for (uint i = 0; i<num_sample; i++)
		{
			// for all samples compute the index in the table based on their node id and class.
			// for all samples count classes in table
			tbl_idx[i] = ((node_id[i] * ord_class) + sample_class[i]);
		}

		// TO BE IMPROVED using SSE instruction
		for (uint i = 0; i<num_sample; i++)
		{
			// for all samples compute the index in the table based on their node id and class.
			// for all samples count classes in table
			table[tbl_idx[i]]++;
		}

		// compute purity of the sample at the current depth of tree
		tree[depth].purity = 0;

		// for each node
		for (uint i = 0; i<num_node_next; i++)
		{
			uint start_idx = i * ord_class;
			uint end_idx = start_idx + ord_class;
			// compute total number of samples in each node
			double sum = 0;
			for (uint j = start_idx; j<end_idx; j++)
				sum += table[j];

			// if there is no sample in the node we skip the node
			if (sum > 0)
			{
				// compute purity of the node
				double node_purity = 0;
				for (uint j = start_idx; j<end_idx; j++)
				{
					double fij = table[j] / sum; // fraction of sample in node i with class j;
					node_purity += (fij * fij);
				}
				// computed weighted average of purity for all nodes
				tree[depth].purity += (sum / num_sample) * node_purity;
			}
		}
		
		PrintfD("\n var_id %u Purity %f", tree[depth].var_id, tree[depth].purity);

		// the number of node is changed now
		num_node = num_node_next;

		return;
	}
};

void UniqueCombination(uint setSize, uint subsetSize, char printTree) 
{
	printf("\n >> %d\t%d", setSize, subsetSize);

	char *table = new char[setSize*setSize];
	memset(table, 0, setSize*setSize);
#define TABLE(I,J)	table[(I*setSize)+J]
	uint *subset = new uint[subsetSize];

	uint treeCNT = 0;
	uint ofs = 0;
	uint depth = 0;

	while (true)
	{
		uint i = ofs;
		for (; i < setSize;)
		{
			// allocate the first element in the subset
			if (depth == 0)
			{
				subset[depth] = i;
				depth++;
				i++;
				continue;
			}

			// for each of the next item first check if repeat exist
			char flag = 0;
			for (uint j = 0; j < depth; j++)
				flag |= TABLE(subset[j], i);
			// if repeat exist check for next element
			if (flag)
			{
				i++;
				continue;
			}
			else // if there where no repeat
			{
				// flag all new repeats
				for (uint j = 0; j < depth; j++)
					TABLE(subset[j], i) = 1;

				// add new element into subset
				subset[depth] = i;
				i++;
				depth++;
				// if enough three is filled print it out
				if (depth == subsetSize)
				{
					treeCNT++;
					if (printTree == 'Y')
					{
						printf("\n");
						for (uint k = 0; k < depth; k++)
							printf("%u\t", subset[k]);
					}
					//and also reset for building new tree
					depth = 0;
					i = ofs;
				}
			}
		}
		// when we exit from this loop it means the first element (ofs) have repeat withe all other element (except afew that can not make a compelte tree)
		// print the last tree if exist
		if (depth > 1)
		{
			treeCNT++;
			if (printTree == 'Y')
			{
				printf("\n");
				for (uint k = 0; k < depth; k++)
					printf("%u\t", subset[k]);
			}
		}
		ofs++;
		depth = 0;
		if (ofs >= setSize)
			break;
	}
	printf("\n>>> Tree Count: %u\n", treeCNT);
	delete[] table;
}

ulint power(ulint b, ulint p)
{
	ulint ret = 1;
	for (ulint i = 0; i < p; i++)
	{
		ret *= b;
	}
	return ret;
}

void ComputeNumberOfPossibleTree(ulint setSize, ulint subsetSize)
{
	ulint i = 0;
	ulint CompeleteTree = 0;
	ulint InCompeleteTree = 0;

	ulint n = setSize;
	ulint m = subsetSize;
	printf("\ni\tA\tterm\t(A*term)\tSub\tTree");
	while (true)
	{
		ulint inSigma = 0;
		for (int j = 0; j <= i; j++)
			inSigma += power(m - 1, j);

		ulint termx = (n > inSigma) ? (n - inSigma) : n;
		ulint term = termx / (m - 1);

		int incompeleteFlag;
		if ((term * (m - 1)) != termx)
			incompeleteFlag = 1;

		ulint A = power(m - 1, i);

		if (n > A)
		{
			CompeleteTree += A * term;
			if (incompeleteFlag)
				InCompeleteTree += A;
			n -= A;
			printf("\n%10llu\t%10llu\t%10llu\t%10llu\t%10llu\t%10llu", i, A, term, (A*term), inSigma, CompeleteTree);
		}
		else
		{
			CompeleteTree += n * term;
			if (incompeleteFlag)
				InCompeleteTree += n;
			printf("\n%10llu\t%10llu\t%10llu\t%10llu\t%10llu\t%10llu", i, n, term, (n*term), inSigma, CompeleteTree);
			break;
		}
		i++;
	}
	printf("\nTotal");
	printf("\n%llu\t%llu\t%llu\t%llu", setSize, subsetSize, CompeleteTree, InCompeleteTree);
}

void PrintHelp(char *progName)
{
	printf(">>> Usage: %s G/A/U/C/P/T\n", progName);
	printf(">>> G: random dataset Generation\n");
	printf(">>> A: Analysis\n");
	printf(">>> U: Unique combinations\n");
	printf(">>> C: Compute Number of possible tree\n");
	printf(">>> P: Compute purity of all variables\n");
	printf(">>> T: Load table (VCF2Table.sh) and simulate uncorrolated variants \n");
	printf(">>> t: Load table (VCF2Table.sh) and simulate corrolated variants \n");
	printf(">>> X: Load table (VCF2Table.sh) and lable \n");
	return;
}

int main(int argc, char *argv[])
{

	if (argc<2)
	{
		PrintHelp(argv[0]);
		return 0;
	}

	srand(time(NULL));

	switch (argv[1][0])
	{
		case 'G':
		{
			if (argc<4)
			{
				printf(">>> Usage: %s G num_var num_samples prefix(ord_var is 3, ord_class is 2)\n", argv[0]);
				return 0;
			}
			Printf("\nSimulation Started.");
			uint num_var = atoi(argv[2]);
			uint num_sample = atoi(argv[3]);
			UNCOMPRESSED_MEMORY vds;
			vds.FillRandom(num_var, num_sample, 3);
			vds.ComputeSamplePurity();
			Printf("\nComputed Sample Purity: %10.9f", vds.sample_purity);
			vds.StoreToFile(argv[4]);
			Printf("\nSimulation End.");
			break;
		}
		case 'T':
		{
			if (argc<5)
			{
				printf(">>> Usage: %s T table_prefix out_prefix num_associated_var (table generated by VCF2Table.sh)\n", argv[0]);
				return 0;
			}
			UNCOMPRESSED_MEMORY vds;
			Printf("\nParse table file.");
			vds.ParseTable(argv[2]);
			Printf("\nSimulation of phenotype.");
			vds.SimLabelNoCorrolation(atoi(argv[4]), argv[3]);
			vds.ComputeSamplePurity();
			Printf("\nComputed Sample Purity: %10.9f", vds.sample_purity);
			vds.StoreToFile(argv[3]);
			Printf("\nSimulation End.");
			break;
		}
		case 't':
		{
			if (argc<5)
			{
				printf(">>> Usage: %s t table_prefix out_prefix num_associated_var (table generated by VCF2Table.sh)\n", argv[0]);
				return 0;
			}
			UNCOMPRESSED_MEMORY vds;
			Printf("\nParse table file.");
			vds.ParseTable(argv[2]);
			Printf("\nSimulation of phenotype.");
			vds.SimLabelCorrolation(atoi(argv[4]), argv[3]);
			vds.ComputeSamplePurity();
			Printf("\nComputed Sample Purity: %10.9f", vds.sample_purity);
			vds.StoreToFile(argv[3]);
			Printf("\nSimulation End.");
			break;
		}
		case 'X':
		{
			if (argc<5)
			{
				printf(">>> Usage: %s X table_prefix out_prefix lable_file_name (table generated by VCF2Table.sh)\n", argv[0]);
				return 0;
			}
			UNCOMPRESSED_MEMORY vds;
			Printf("\nParse table file.");
			vds.ParseTable(argv[2]);
			Printf("\nLoading phenotype.");
			vds.Loadlabel(argv[4]);
			vds.ComputeSamplePurity();
			Printf("\nComputed Sample Purity: %10.9f", vds.sample_purity);
			vds.StoreToFile(argv[3]);
			Printf("\nSimulation End.");
			break;
		}
		case 'A':
		{
			if (argc<7)
			{
				printf(">>> Usage: %s A prefix max_depth tree_per_file num_tree tree_file_prefix\n", argv[0]);
				return 0;
			}
			UNCOMPRESSED_MEMORY vds;
			vds.LoadFromFile(argv[2]);
			CORE core;
			uint thread_id = 0;
			core.Init(vds.num_sample, vds.sample_class, 2, vds.sample_purity, atoi(argv[3]), atoi(argv[4]), vds.num_var, 3, memory, &vds, argv[6], thread_id);
			Printf("\nBuilding Forest ");
			core.BuildForest(atoi(argv[5]));
			break;
		}
		case 'U':
		{
			if (argc<5)
			{
				printf(">>> Usage: %s U setSize subsetSize PrintTree(Y/N)\n", argv[0]);
				return 0;
			}
			UniqueCombination(atoi(argv[2]), atoi(argv[3]), argv[4][0]);
			break;
		}
		case 'C':
		{
			if (argc<4)
			{
				printf(">>> Usage: %s C setSize subsetSize\n", argv[0]);
				return 0;
			}
			ComputeNumberOfPossibleTree(atoi(argv[2]), atoi(argv[3]));
			break;
		}
		case 'P':
		{
			if (argc<3)
			{
				printf(">>> Usage: %s P file_name max_depth tree_per_file num_tree tree_file_prefix\n", argv[0]);
				return 0;
			}
			UNCOMPRESSED_MEMORY vds;
			vds.LoadFromFile(argv[2]);
			CORE core;
			uint thread_id = 0;
			char dummy[10];
			strcpy(dummy, "DUMMY");
			core.Init(vds.num_sample, vds.sample_class, 2, vds.sample_purity, 3, 0, vds.num_var, 3, memory, &vds, dummy, thread_id);
			Printf("\nComputing purity of all variable");
			core.ComputeAndPrintAllInfoGained();
			break;
		}
		default:
		{
			PrintHelp(argv[0]);
		}
	}
	Printf("\n");
	return 0;
}