// ==========================================================================
//                                 build_filter.cpp
// ==========================================================================
// Copyright (c) 2017-2022, Temesgen H. Dadi, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Temesgen H. Dadi or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL TEMESGEN H. DADI OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Temesgen H. Dadi <temesgen.dadi@fu-berlin.de>
// ==========================================================================

#define BUILD_FILTER
// ----------------------------------------------------------------------------
// STL headers
// ----------------------------------------------------------------------------
#include <string>
#include <vector>
#include <mutex>
#include <condition_variable>
#include <future>
#include <thread>

// ----------------------------------------------------------------------------
// SeqAn headers
// ----------------------------------------------------------------------------
#include <seqan/index.h>

// ----------------------------------------------------------------------------
// App headers
// ----------------------------------------------------------------------------
#include "store_seqs.h"
#include "misc_timer.h"
#include "misc_types.h"
#include "bits_matches.h"
#include "misc_options.h"
#include "misc_options_dis.h"
#include "kdx_filter.h"
#include "bloom_filter.h"
#include "index_fm.h"


using namespace seqan;

// ----------------------------------------------------------------------------
// Class Options
// ----------------------------------------------------------------------------

struct Options
{
    CharString      kmersDir;
    CharString      bigIndexDir;
    CharString      indicesDir;
    CharString      uniqeKmerDir;
    CharString      filterFile;

    uint32_t        kmerSize;
    uint32_t        numberOfBins;
    uint64_t        bloomFilterSize;
    uint32_t        numberOfHashes;
    unsigned        threadsCount;
    bool            verbose;

    FilterType      filterType;

    std::vector<std::string> filterTypeList;

    Options() :
    kmerSize(20),
    numberOfBins(64),
    bloomFilterSize(8589934592 + filterMetadataSize), // 1GB
    numberOfHashes(4),
    threadsCount(1),
    verbose(false),
    filterType(BLOOM),
    filterTypeList({"bloom", "kmer_direct", "none"})
    {
    }
};

// ==========================================================================
// Functions
// ==========================================================================

// ----------------------------------------------------------------------------
// Function setupArgumentParser()
// ----------------------------------------------------------------------------

void setupArgumentParser(ArgumentParser & parser, Options const & options)
{
    setAppName(parser, "yara_build_filter");
    setShortDescription(parser, "Build Filter for Distributed Yara");
    setCategory(parser, "Read Mapping");

    setDateAndVersion(parser);
    setDescription(parser);

    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIREFERENCE FILES DIR \\fP>");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_PREFIX, "REFERENCE FILE DIR"));
    //    setValidValues(parser, 0, SeqFileIn::getFileExtensions());
    setHelpText(parser, 0, "A directory containing reference genome files.");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_PREFIX, "INDICES DIR"));
    setHelpText(parser, 1, "A directory containing fm indices.");

    addArgument(parser, ArgParseArgument(ArgParseArgument::OUTPUT_PREFIX, "UNIQUE KMERS"));
    setHelpText(parser, 2, "A directory to save unique kmers to.");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_PREFIX, "BIG INDICES DIR"));
    setHelpText(parser, 3, "A directory containing the big indices");

    addOption(parser, ArgParseOption("v", "verbose", "Displays verbose output."));

    addSection(parser, "Output Options");

    addOption(parser, ArgParseOption("o", "output-file", "Specify an output filename for the filter. \
                                     Default: use the directory name of reference genomes.", ArgParseOption::OUTPUT_FILE));
    setValidValues(parser, "output-file", "filter");

    addOption(parser, ArgParseOption("b", "number-of-bins", "The number of bins (indices) for distributed mapper",
                                     ArgParseOption::INTEGER));

    setMinValue(parser, "number-of-bins", "1");
    setMaxValue(parser, "number-of-bins", "10000");

    addOption(parser, ArgParseOption("t", "threads", "Specify the number of threads to use.", ArgParseOption::INTEGER));
    setMinValue(parser, "threads", "1");
    setMaxValue(parser, "threads", "2048");
    setDefaultValue(parser, "threads", options.threadsCount);

    addOption(parser, ArgParseOption("ft", "filter-type", "type of filter to build",
                                     ArgParseOption::STRING));
    setValidValues(parser, "filter-type", options.filterTypeList);
    setDefaultValue(parser, "filter-type",  options.filterTypeList[options.filterType]);


    addOption(parser, ArgParseOption("k", "kmer-size", "The size of kmers for bloom_filter",
                                     ArgParseOption::INTEGER));
    setMinValue(parser, "kmer-size", "14");
    setMaxValue(parser, "kmer-size", "32");

    addOption(parser, ArgParseOption("nh", "num-hash", "Specify the number of hash functions to use for the bloom filter.", ArgParseOption::INTEGER));
    setMinValue(parser, "num-hash", "2");
    setMaxValue(parser, "num-hash", "5");
    setDefaultValue(parser, "num-hash", options.numberOfHashes);

    addOption(parser, ArgParseOption("bs", "bloom-size", "The size of bloom filter in GB.", ArgParseOption::INTEGER));
    setMinValue(parser, "bloom-size", "1");
    setMaxValue(parser, "bloom-size", "512");
    setDefaultValue(parser, "bloom-size", 1);
}

// ----------------------------------------------------------------------------
// Function parseCommandLine()
// ----------------------------------------------------------------------------

ArgumentParser::ParseResult
parseCommandLine(Options & options, ArgumentParser & parser, int argc, char const ** argv)
{
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Parse verbose output option.
    getOptionValue(options.verbose, parser, "verbose");

    // Parse contigs input file.
    getArgumentValue(options.kmersDir, parser, 0);
    getArgumentValue(options.indicesDir, parser, 1);
    getArgumentValue(options.uniqeKmerDir, parser, 2);
    getArgumentValue(options.bigIndexDir, parser, 3);

    // Parse contigs index prefix.
    getOptionValue(options.filterFile, parser, "output-file");
    if (!isSet(parser, "output-file"))
    {
        options.filterFile = trimExtension(options.kmersDir);
        append(options.filterFile, "bloom.filter");
    }

    getOptionValue(options.filterType, parser, "filter-type", options.filterTypeList);

    if (isSet(parser, "number-of-bins")) getOptionValue(options.numberOfBins, parser, "number-of-bins");
    if (isSet(parser, "kmer-size")) getOptionValue(options.kmerSize, parser, "kmer-size");
    if (isSet(parser, "threads")) getOptionValue(options.threadsCount, parser, "threads");
    if (isSet(parser, "num-hash")) getOptionValue(options.numberOfHashes, parser, "num-hash");

    uint64_t bloomSize;
    if (getOptionValue(bloomSize, parser, "bloom-size"))
    {
        if ((bloomSize & (bloomSize - 1)) == 0)
        {
            options.bloomFilterSize = bloomSize * 8589934592 + filterMetadataSize; // 8589934592 = 1GB
        }
        else
        {
            std::cerr <<"[ERROR] --bloom-size (-bs) parameter should be a power of 2!" << std::endl;
            exit(1);
        }
    }
    return ArgumentParser::PARSE_OK;
}


// ----------------------------------------------------------------------------
// Function get_unique_kmers()
// ----------------------------------------------------------------------------
inline void get_unique_kmers(Options & options)
{
    typedef YaraFMConfig<uint16_t, uint32_t, uint64_t>    TIndexConfig;
//    typedef YaraFMConfig<uint8_t, uint16_t, uint32_t>    TIndexConfig;
    typedef FMIndex<void, TIndexConfig>                             TIndexSpec;
    typedef Index<typename TIndexConfig::Text, TIndexSpec>          TIndex;


    String<uint64_t> limits;

    CharString contigsLimitFile(options.bigIndexDir);
    append(contigsLimitFile, ".txt.size");
    open(limits, toCString(contigsLimitFile), OPEN_RDONLY);

    std::cout << limits[1] << ", "  << limits[0] << ", "  << limits[2] << "\n";
    std::string comExt = commonExtension(options.kmersDir, options.numberOfBins);
    std::map<CharString, uint32_t> bin_map;
    std::map<uint32_t, uint32_t> big_bin_map;

    std::vector<uint64_t> uniq_counts(options.numberOfBins, 0);
    std::vector<uint64_t> kmer_counts(options.numberOfBins, 0);

    typedef SeqStore<void, YaraContigsConfig<> >                    TContigs;

    for (uint32_t binNo = 0; binNo < options.numberOfBins; ++binNo)
    {
        CharString fm_index_file;
        appendFileName(fm_index_file, options.indicesDir, binNo);

        TContigs tmpContigs;

        if (!open(tmpContigs, toCString(fm_index_file), OPEN_RDONLY))
            throw RuntimeError("Error while opening reference file.");

        for (uint32_t i = 0; i < length(tmpContigs.names); ++i)
        {
            CharString s = (CharString)tmpContigs.names[i];
            bin_map[s] = binNo;
//            std::cout << "small i = " << i  << "Bin " << binNo << "\tRef Name " << s << std::endl ;
        }
    }

    TContigs allContigs;

    if (!open(allContigs, toCString(options.bigIndexDir), OPEN_RDONLY))
        throw RuntimeError("Error while opening reference file.");

    for (uint32_t i = 0; i < length(allContigs.names); ++i)
    {
        CharString s = (CharString)allContigs.names[i];
        big_bin_map[i] = bin_map[s];
//        if(bin_map[s] == 0)
//        std::cout << "Big i = " << i  << "Bin " << bin_map[s] << "\tRef Name " << s << std::endl ;
    }

    TIndex big_fm_index;
    if (!open(big_fm_index, toCString(options.bigIndexDir), OPEN_RDONLY))
        throw "ERROR: Could not open the index.";


    Semaphore thread_limiter(options.threadsCount);
    std::vector<std::future<void>> tasks;

    for (uint32_t binNo = 0; binNo < options.numberOfBins; ++binNo)
    {
        tasks.emplace_back(std::async([=, &thread_limiter, &uniq_counts, &kmer_counts, &big_fm_index] {
        Critical_section _(thread_limiter);

        uint32_t batchSize = 10000;
        Finder<TIndex> finder(big_fm_index);

//        CharString rawUniqKmerFile;
//        appendFileName(rawUniqKmerFile, options.uniqeKmerDir, binNo);
//        append(rawUniqKmerFile, ".fa.gz");
//
//        SeqFileOut rawFileOut(toCString(rawUniqKmerFile));


        CharString fastaFile;
        appendFileName(fastaFile, options.kmersDir, binNo);
        append(fastaFile, comExt);
        SeqFileIn seqFileIn;
        if (!open(seqFileIn, toCString(fastaFile)))
        {
            std::cerr <<"Unable to open contigs File: " << toCString(fastaFile) << std::endl;
            exit(1);
        }
        StringSet<CharString> ids;
        StringSet<IupacString> seqs;

        while(!atEnd(seqFileIn))
        {
//            CharString id;
//            IupacString seq;
//            readRecord(id, seq, seqFileIn);
//            ++counter;
            readRecords(ids, seqs, seqFileIn, batchSize);
            uint32_t len = length(seqs);
            kmer_counts[binNo] += len;
            for(uint32_t i = 0; i<len; ++i)
            {
                bool uniq = true;
                while (find(finder, seqs[i]))
                {
                    auto pos = position(finder);
                    uint32_t rID = getValueI1(pos);
                    auto bi = big_bin_map.find(rID);

                    if(bi != big_bin_map.end() && binNo != bi->second)
                    {
                        uniq = false;
                        break;
                    }
                }
                if (uniq)
                {
//                    writeRecord(rawFileOut, counter, seqs[i]);
                    ++uniq_counts[binNo];
                }
                goBegin(finder);    // move Finder to the beginning of the text
                clear(finder);      // reset Finder
            }
            clear(ids);
            clear(seqs);
        }
        close(seqFileIn);
//        close(rawFileOut);
//        std::cerr << counter <<" kmers from bin " << binNo << std::endl;
//        std::cerr << uniq_counter <<" unique kmers from bin " << binNo << std::endl;
        }));
    }
    for (auto &&task : tasks)
    {
        task.get();
    }

    std::cout << "binNo\tkmer count\tuniq count\n";
    for (uint32_t binNo = 0; binNo < options.numberOfBins; ++binNo)
    {
        std::cout << binNo  << "\t" <<  kmer_counts[binNo] << "\t" << uniq_counts[binNo] << std::endl ;
    }
}

// ----------------------------------------------------------------------------
// Function build_filter()
// ----------------------------------------------------------------------------
template <typename TFilter>
inline void build_filter(Options & options, TFilter & filter)
{
    std::string comExt = commonExtension(options.kmersDir, options.numberOfBins);

    Semaphore thread_limiter(options.threadsCount);
    std::vector<std::future<void>> tasks;

    Timer<double>       timer;
    Timer<double>       globalTimer;
    start (timer);
    start (globalTimer);

    for (uint32_t binNo = 0; binNo < options.numberOfBins; ++binNo)
    {
        tasks.emplace_back(std::async([=, &thread_limiter, &filter] {
            Critical_section _(thread_limiter);

            Timer<double>       binTimer;
            start (binTimer);

            CharString fastaFile;
            appendFileName(fastaFile, options.kmersDir, binNo);
            append(fastaFile, comExt);

            filter.addFastaFile(fastaFile, binNo);

            stop(binTimer);
            if (options.verbose)
            {
                mtx.lock();
                std::cerr <<"[bin " << binNo << "] Done adding kmers!\t\t\t" << binTimer << std::endl;
                mtx.unlock();
            }
        }));
    }
    for (auto &&task : tasks)
    {
        task.get();
    }
    stop(timer);
    if (options.verbose)
        std::cerr <<"All bins are done adding kmers!\t\t" << timer << std::endl;

    start(timer);
    filter.save(toCString(options.filterFile));
    stop(timer);
    if (options.verbose)
        std::cerr <<"Done saving filter (" << filter.size_mb() <<" MB)\t\t" << timer << std::endl;
    stop(globalTimer);
    std::cerr <<"\nFinshed in \t\t\t" << globalTimer << std::endl;

}


// ----------------------------------------------------------------------------
// Function main()
// ----------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
    ArgumentParser parser;
    Options options;
    setupArgumentParser(parser, options);

    ArgumentParser::ParseResult res = parseCommandLine(options, parser, argc, argv);

    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    try
    {
        get_unique_kmers(options);
//        if (options.filterType == BLOOM)
//        {
//            SeqAnBloomFilter<> filter  (options.numberOfBins,
//                                        options.numberOfHashes,
//                                        options.kmerSize,
//                                        options.bloomFilterSize);
//
//            build_filter(options, filter);
//        }
//        else if (options.filterType == KMER_DIRECT)
//        {
//            uint64_t vec_size = (1u << (2 * options.kmerSize));
//            vec_size *= options.numberOfBins;
//            vec_size += filterMetadataSize;
//            SeqAnKDXFilter<> filter (options.numberOfBins, options.kmerSize, vec_size);
//
//            build_filter(options, filter);
//        }

    }
    catch (Exception const & e)
    {
        std::cerr << getAppName(parser) << ": " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
