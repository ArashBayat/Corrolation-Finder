/*
* Core.cpp
*
*  Created on: 1 Feb. 2018
*      Author: bay041
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

typedef enum
{
	memory,					// variable data stored in a table in the memory
	compressed_memory,		// variable data table is compressed in the memory
	RASM						// random access to records with sequential access to file
} MEMORY_PLAN;

typedef struct T_UNCOMPRESSED_MEMORY
{
	SAMPLE_COUNTER num_sample;	// number of samples
	uint num_var;				// number of variables

	VAR_VALUE *data;		// keep value for each variable
	VAR_ORD *ord;		// ordinallity of each variable
	SAMPLE_CLASS *sample_class; // class labels for samples

	VAR_VALUE *Get_Var_Value(uint index)
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

		// fill random value for variables
		for (uint i = 0; i<num_var; i++)
		{
			ord[i] = a_ord;
			for (uint j = 0; j<num_sample; j++)
				data[(i * num_sample) + j] = rand() % ord[i];
		}

		// assign random class for samples
		for (uint i = 0; i<num_sample; i++)
			sample_class[i] = rand() % 2;
	}

	void StoreToFile(const char *file_name)
	{
		FILE *file = fopen(file_name, "wb");
		NULL_CHECK(file);
		fwrite(this, sizeof(T_UNCOMPRESSED_MEMORY), 1, file);
		fwrite(this->data, sizeof(VAR_VALUE), (num_var * num_sample), file);
		fwrite(this->ord, sizeof(VAR_ORD), (num_var), file);
		fwrite(this->sample_class, sizeof(SAMPLE_CLASS), (num_sample), file);
	}

	void LoadFromFile(const char *file_name)
	{
		FILE *file = fopen(file_name, "rb");
		NULL_CHECK(file);
		fread(this, sizeof(T_UNCOMPRESSED_MEMORY), 1, file);

		data = new VAR_VALUE[num_var * num_sample];
		NULL_CHECK(data);
		ord = new VAR_ORD[num_var];
		NULL_CHECK(ord);
		sample_class = new SAMPLE_CLASS[num_sample];
		NULL_CHECK(sample_class);

		fread(this->data, sizeof(VAR_VALUE), (num_var * num_sample), file);
		fread(this->ord, sizeof(VAR_ORD), (num_var), file);
		fread(this->sample_class, sizeof(SAMPLE_CLASS), (num_sample), file);
	}

	~T_UNCOMPRESSED_MEMORY()
	{
		delete[] data;
		delete[] ord;
		delete[] sample_class;
	}
} UNCOMPRESSED_MEMORY;

typedef struct
{
	uint ord_var;	// tree[i].ord_var=j indicates that variable used in depth 'i' to split data has the ordinality of 'j'. Thus the value of the variable is between 0 and (ord_var[i]-1)
	uint var_id;		// tree[i].var_id=j indicates that id of variable used in depth 'i' to split samples is 'j'.
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

		if (file_tree) // not to leave output file open
			fclose(file_tree);
	}

	void Init(SAMPLE_COUNTER a_num_sample, SAMPLE_CLASS *a_sample_class, uint a_ord_class, uint a_max_depth, uint a_num_tree_per_file, uint a_num_var, uint a_max_var_ord, MEMORY_PLAN a_memory_plan, void *a_dataset, char *prefix, uint a_thread_id)
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

		if (strlen(prefix) >= 4096)
			Error("Tree prefix file name is too long. Use shorter path (max 4095 characters)");
		strcpy(tree_file_name_prefix, prefix);
		thread_id = a_thread_id;

		num_var = a_num_var;
		memory_plan = a_memory_plan;
		switch (memory_plan)
		{
		case memory:
			uncompressed_memory = (UNCOMPRESSED_MEMORY *)a_dataset;
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
		for (int i = 0; i<num_tree; i++)
			BuildTree();
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

	void TakeRandomVariable()
	{
		// select a random variable
		uint rnd = (uint)(rand());
		rnd *= (uint)(rand());
		rnd %= num_var;

		// based on memory plan load variable data
		switch (memory_plan)
		{
		case memory:
			tree[depth].var_id = rnd;
			tree[depth].ord_var = uncompressed_memory->ord[rnd];
			var_value = uncompressed_memory->Get_Var_Value(rnd);
			break;
		case compressed_memory:
			break;
		case RASM:
			break;
		}

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
		}
		// finally store the tree
		StoreTree();
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
//#define NO_TBALE
#ifndef NO_TABLE
		// TO BE IMPROVED using SSE instruction
		for (uint i = 0; i<num_sample; i++)
		{
			// for all samples compute the index in the table based on their node id and class.
			// for all samples count classes in table
			table[tbl_idx[i]]++;
		}
#endif
		// the number of node is changed now
		num_node = num_node_next;

		// compute purity of the sample at the current depth of tree
		tree[depth].purity = 0;
#ifdef NO_TABLE
		for (uint i = 0; i < num_node; i++)
			tree[depth].purity += node_id[i];
#endif
		// for each node
#ifndef NO_TABLE
		for (uint i = 0; i<num_node; i++)
		{
			uint start_idx = i * ord_class;
			uint end_idx = start_idx + ord_class;
			// compute total number of samples in each node
			double sum = 0;
			for (int j = start_idx; j<end_idx; j++)
				sum += table[j];

			// if there is no sample in the node we skip the node
			if (sum > 0)
			{
				// compute purity of the node
				double node_purity = 0;
				for (int j = start_idx; j<end_idx; j++)
				{
					double fij = table[j] / sum; // fraction of sample in node i with class j;
					node_purity += (fij * fij);
				}
				// computed weighted average of purity for all nodes
				tree[depth].purity += (sum / num_sample) * node_purity;
			}
		}
#endif
		PrintfD("\n var_id %u Purity %f", tree[depth].var_id, tree[depth].purity);
		// increase the depth of the tree
		depth++;

		return;
	}
};

int main(int argc, char *argv[])
{
	if (argc<2)
	{
		printf(">>> Usage: %s R/A (R for random dataset Generation/ A for analysis)\n", argv[0]);
		return 0;
	}
	switch (argv[1][0])
	{
	case 'G':
	{
		if (argc<4)
		{
			printf(">>> Usage: %s R num_var num_samples file_name(ord_var is 3, ord_class is 2)\n", argv[0]);
			return 0;
		}
		Printf("\nSimulation Started.");
		uint num_var = atoi(argv[2]);
		uint num_sample = atoi(argv[3]);
		UNCOMPRESSED_MEMORY vds;
		vds.FillRandom(num_var, num_sample, 3);
		vds.StoreToFile(argv[4]);
		Printf("\nSimulation End.");
		break;
	}
	case 'A':
	{
		if (argc<7)
		{
			printf(">>> Usage: %s A file_name max_depth tree_per_file num_tree tree_file_prefix\n", argv[0]);
			return 0;
		}
		UNCOMPRESSED_MEMORY vds;
		vds.LoadFromFile(argv[2]);
		CORE core;
		uint thread_id = 0;
		core.Init(vds.num_sample, vds.sample_class, 2, atoi(argv[3]), atoi(argv[4]), vds.num_var, 3, memory, &vds, argv[6], thread_id);
		Printf("\nBuilding Forest ");
		core.BuildForest(atoi(argv[5]));
		break;
	}
	default:
	{
		printf(">>> Usage: %s R/A (R for random dataset Generation/ A for analysis)\n", argv[0]);
	}
	}
	Printf("\n");
	return 0;
}