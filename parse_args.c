#ifndef PARSE_ARGS_INCLUDED
#define PARSE_ARGS_INCLUDED

#include <argp.h>
#include <stdlib.h>

const char *argp_program_version = "network embedder v0.01";
const char *argp_program_bug_address = "<saska.donges@helsinki.fi>";
static char doc[] = "Simple 2d embedding generator for netwoks.\nCalling should be prefixed with \"mpirun -np ntasks\",\nwhere ntasks is the number of processes to run.";
static char args_doc[] = "FILE";
static struct argp_option options[] = { 
    { 0, 'n', "N", 0, "Highest node id in input. Required"},
    { 0, 'k', "K", 0, "Number of closest and furthest neighbours to consider. Default or 0 to automatically set."},
    { "reverse_penalty", 'p', "P", 0, "Penalty for following reverse edges. Defaults to 3."},
    { "xsize", 'x', "X", 0, "Number of rows in universe. Default or 0 for automatically calculating."},
    { "ysize", 'y', "Y", 0, "Number of columns in universe. Default or 0 for automatically calculating."},
    { "iters", 'i', "I", 0, "Maximum number of iterations to run."},
    { "interval", 't', "T", 0, "Output universe state every T iterations. Default or 0 for no intermediate output."},
    { 0 } 
};

struct arguments {
    int n, k, p, x, y, id, t;
    int iters;
    int interval;
    char *file_path;
};

static void invalid(char *name, int val, struct argp_state *state, int id) {
    if (id == 0) {
        fprintf(stderr, "Invalid value for %s: %d\n", name, val);
        argp_usage(state);
    }
    exit(argp_err_exit_status);
}

static void too_many(struct argp_state *state, int id) {
    if (id == 0) {
        fprintf(stderr, "Too many arguments passed to executable.\n");
        argp_usage(state);
    }
    exit(argp_err_exit_status);
}

static void missing_file(struct argp_state *state, int id)  {
    if (id == 0) {
        fprintf(stderr, "Input file is required.\n");
        argp_usage(state);
    }
    exit(argp_err_exit_status);
}

static error_t parse_opt(int key, char *arg, struct argp_state *state) {
    struct arguments *arguments = state->input;
    switch (key) {
        case ARGP_KEY_ARG:
            if (state->arg_num > 0) too_many(state, arguments->id);
            arguments->file_path = arg;
#ifdef DEBUG
            fprintf(stderr, "Parsed input file path: %s\n", arguments->file_path);
#endif 
            break;
        case ARGP_KEY_END:
            if (state->arg_num != 1) missing_file(state, arguments->id);
            if (arguments->n < 1) invalid("n", arguments->n, state, arguments->id);
            if (arguments->k < 0) invalid("k", arguments->k, state, arguments->id);
            if (arguments->x < 0) invalid("x", arguments->x, state, arguments->id);
            if (arguments->y < 0) invalid("y", arguments->y, state, arguments->id);
            if (arguments->p < 0) invalid("p", arguments->p, state, arguments->id);
            if (arguments->iters < 1) invalid("i", arguments->iters, state, arguments->id);
            if (arguments->interval < 0) invalid("t", arguments->interval, state, arguments->id);
            break;
        case 'n': 
            arguments->n = atoi(arg);
#ifdef DEBUG
            fprintf(stderr, "Parsed n = %d\n", arguments->n);
#endif 
            break;
        case 'k': 
            arguments->k = atoi(arg); 
#ifdef DEBUG
            fprintf(stderr, "Parsed k = %d\n", arguments->k);
#endif 
            break;
        case 'p': 
            arguments->p = atoi(arg);
#ifdef DEBUG
            fprintf(stderr, "Parsed p = %d\n", arguments->p);
#endif 
            break;
        case 'x':
            arguments->x = atoi(arg);
#ifdef DEBUG
            fprintf(stderr, "Parsed x = %d\n", arguments->x);
#endif 
            break;
        case 'y':
            arguments->y = atoi(arg);
#ifdef DEBUG
            fprintf(stderr, "Parsed y = %d\n", arguments->y);
#endif 
            break;
        case 'i':
            arguments->iters = atoi(arg);
#ifdef DEBUG
            fprintf(stderr, "Parsed i = %d\n", arguments->iters);
#endif 
            break;
        case 't':
            arguments->interval = atoi(arg);
#ifdef DEBUG
            fprintf(stderr, "Parsed t = %d\n", arguments->interval);
#endif 
            break;
        default: 
            return ARGP_ERR_UNKNOWN;
    }   
    return 0;
}

static struct argp argp = { options, parse_opt, args_doc, doc };

struct arguments *parse_arguments(int argc, char *argv[], int id) {
    struct arguments *args = calloc(1, sizeof(struct arguments));
    args->id = id;
    args->iters = 1000;
    argp_parse(&argp, argc, argv, 0, 0, args);
    return args;
}

#endif