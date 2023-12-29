#pragma once

#ifdef TRACES

#include <chrono>
#include <iostream>
#include <fstream>
#include <vector>

#ifdef OMP_USE

#include <omp.h>

#endif

#define CSV_HEADER "Gen,Resource,Stamp,TimeStamp_Start,TimeStamp_End,Duration,Indiv_id"
#define GET_TIME std::chrono::steady_clock::now().time_since_epoch().count()

#ifdef OMP_USE
#define GET_RESOURCE omp_get_thread_num()
#define GET_MAX_RESOURCES omp_get_max_threads()
#else
#define GET_RESOURCE 0
#define GET_MAX_RESOURCES 1
#endif

namespace time_tracer {
    std::vector<const char *> stamp_name;

    static std::ofstream trace_file;
    static std::vector<long> *starts = nullptr;
    static std::vector<long> *ends = nullptr;
    static std::vector<int> *stamp_int = nullptr;
    static std::vector<int> *indiv_id = nullptr;

    static int nb_resources = 0;

    static void init_tracer(const char *trace_file_name, std::vector<const char *> stamps) {
        trace_file.open(trace_file_name, std::ofstream::trunc);
        trace_file << CSV_HEADER << std::endl;

        stamp_name = std::move(stamps);
        nb_resources = GET_MAX_RESOURCES;

        starts = new std::vector<long>[nb_resources];
        ends = new std::vector<long>[nb_resources];
        stamp_int = new std::vector<int>[nb_resources];
        indiv_id = new std::vector<int>[nb_resources];
    }

    static void stop_tracer() {
        delete[] starts;
        delete[] ends;
        delete[] stamp_int;
        starts = nullptr;
        ends = nullptr;
        stamp_int = nullptr;
        indiv_id = nullptr;
        nb_resources = 0;
        stamp_name.clear();

        trace_file.close();
    }

    static void timestamp_start() {
        starts[GET_RESOURCE].push_back(GET_TIME);
    }

    static void timestamp_end(int stamp_id) {
        auto resource = GET_RESOURCE;
        ends[resource].push_back(GET_TIME);
        stamp_int[resource].push_back(stamp_id);
        indiv_id[resource].push_back(-1);
    }

    static void timestamp_end(int stamp_id, int id) {
        auto resource = GET_RESOURCE;
        ends[resource].push_back(GET_TIME);
        stamp_int[resource].push_back(stamp_id);
        indiv_id[resource].push_back(id);
    }

    static void set_traces() {
        for (int i = 0; i < nb_resources; ++i) {
            starts[i].clear();
            ends[i].clear();
            stamp_int[i].clear();
        }
    }

    static void write_traces(int generation) {
        for (int res = 0; res < nb_resources; ++res) {
            for (int stamp = 0; stamp < stamp_int[res].size(); ++stamp) {
                trace_file << generation
                           << "," << (res + 1)
                           << "," << stamp_name[stamp_int[res][stamp]]
                           << "," << starts[res][stamp]
                           << "," << ends[res][stamp]
                           << "," << ends[res][stamp] - starts[res][stamp]
                           << "," << indiv_id[res][stamp]
                           << std::endl;
            }
        }
        trace_file.flush();
        set_traces();
    }
}

#define INIT_TRACER(file_name, ...) time_tracer::init_tracer(file_name, __VA_ARGS__);

#define TIMESTAMP(STAMP, BLOCK) { \
time_tracer::timestamp_start(); \
BLOCK \
time_tracer::timestamp_end(STAMP); \
}

#define TIMESTAMP_ID(STAMP, ID, BLOCK) { \
time_tracer::timestamp_start(); \
BLOCK \
time_tracer::timestamp_end(STAMP, ID); \
}

#define FLUSH_TRACES(generation) time_tracer::write_traces(generation);
#define STOP_TRACER time_tracer::stop_tracer();

#else //#ifndef TREACES

#define INIT_TRACER(file_name, ...)
#define TIMESTAMP(STAMP, BLOCK) BLOCK
#define TIMESTAMP_ID(STAMP, ID, BLOCK) BLOCK
#define FLUSH_TRACES(generation)
#define STOP_TRACER

#endif //TRACES