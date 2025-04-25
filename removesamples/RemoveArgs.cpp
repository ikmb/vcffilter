/*
 *    Copyright (C) 2025 by Lars Wienbrandt,
 *    Institute of Clinical Molecular Biology, Kiel University
 *
 *    This file is part of Vcffilter.
 *
 *    Vcffilter is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    Vcffilter is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with Vcffilter. If not, see <https://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <iomanip>
#include <string>

#include <boost/program_options.hpp>

#include "RemoveArgs.h"

namespace bpo = boost::program_options;

using namespace std;
using namespace bpo;

/**
 * @brief Constructs, parses and verifies all command-line arguments
 *
 * This function constructs the Args object and parses all given command-line
 * arguments. If unknown arguments and/or missing obligatory arguments are
 * detected, this function does not return and instead prints an appropriate
 * help message and calls exit(), executing all exit handlers registered
 * up to the parseArgs call.
 *
 * @param argc argument count
 * @param argv argument vector
 * @return an rvalue reference of a new Args object containing all defined or defaulted CLI arguments
 */
/*static*/ RemoveArgs RemoveArgs::parseArgs(int argc, char *argv[]) {
    RemoveArgs args {argc, argv};

    if (args.count("help")) {
        args.printHelp(argv[0], cerr);
        exit(EXIT_SUCCESS);
    }

    if (args.count("version")) {
        args.printVersion(cerr);
        exit(EXIT_SUCCESS);
    }

    if(!args.count("skipidfile")) {
        cerr << "You need to specify a file with sample IDs that should be excluded from the input VCF stream." << endl;
        exit(EXIT_FAILURE);
    }

    // set variables
    args.parseVars();

    // mac and maf filter cannot be used together
    if (args.macfilter && args.maffilter > 0) {
        cerr << "ERROR: MAC and MAF filter cannot be used together." << endl;
        exit(EXIT_FAILURE);
    }

    return args;
}

ostream &operator<<(ostream &out, const RemoveArgs &args) {
    variables_map::const_iterator it;

    long name_width = 0;
    long value_width = 0;
    long orig_width = out.width();

    // collect field widths
    for(it = args.vars.begin(); it != args.vars.end(); ++it) {
        long this_width = static_cast<long>(it->first.length());

        if(this_width > name_width) {
            name_width = this_width;
        }

        this_width = 0;

        if(it->second.value().type() == typeid(string)) {
            this_width = static_cast<long>(it->second.as<string>().length());
        }

        if(it->second.value().type() == typeid(int)) {
            this_width = static_cast<long>(log10(it->second.as<int>()));
        }

        if(this_width > value_width) {
            value_width = this_width;
        }
    }

    // dump option values
    out.setf(ios::adjustfield, ios::left);

    for(it = args.vars.begin(); it != args.vars.end(); ++it) {
        out.width(name_width + 2); // more space for the colon
        out << (it->first + ":") << endl;

        out.width(value_width);

        if(it->second.value().type() == typeid(string)) {
            out << it->second.as<string>() << endl;
        } else if(it->second.value().type() == typeid(int)) {
            out << it->second.as<int>() << endl;
        } else {
            out << "(unknown)" << endl;
        }
    }

    resetiosflags(ios::adjustfield);
    out.width(orig_width);

    return out;
}

RemoveArgs::RemoveArgs(int argc, char *argv[]) :
    opts_regular("Program options"),
    opts_hidden("Hidden options (only visible in debug mode)")
    {

    opts_regular.add_options()
    ("help,h", "produces this help message and exits")
    ("version,v", "prints version information and exits")
    ("macfilter", value<size_t>(&macfilter)->default_value(0), "only variants with a minor allele count >= value are returned")
    ("maffilter", value<float>(&maffilter)->default_value(0.0), "only variants with a minor allele frequency >= value are returned")
    ("missfilter", value<float>(&missfilter)->default_value(0.0), "only variants with a genotype missingness rate < value are returned")
    ;

    opts_hidden.add_options()
    ("skipidfile", value<string>(&skipidfilename), "file with sample IDs that should be removed from the VCF input stream (optional for positional arg #1)")
    ("debug", "produce lots of debug output")
    ;

    opts_positional.add("skipidfile", 1);

    parse(argc, argv);
}

void RemoveArgs::parse(int argc, char *argv[]) {
    bpo::options_description all_options;

    // combine all options
    all_options.add(opts_regular);
    all_options.add(opts_hidden);

    // do the actual parsing
    store(command_line_parser(argc, argv).options(all_options).positional(opts_positional).run(), vars);
    notify(vars);

}

void RemoveArgs::parseVars() {

    if (vars.count("debug"))
        debug = true;

}

bool RemoveArgs::isDefined(const string &optname) const {
    bool found = false;
    found = !!this->opts_regular.find_nothrow(optname, false); // return null (-> false) if option has not been found
    found |= !!this->opts_hidden.find_nothrow(optname, false);
    return found;
}

void RemoveArgs::printHelp(const string &progname, ostream &out) const {
    out << "Usage: " << progname << " <skipidfile> [options]" << endl << endl;
    out << opts_regular << endl;
    if (debug) {
        out << opts_hidden << endl;
    }
    out << endl;

    out << " The tool reads a VCF file stream from stdin and produces a VCF file stream without the samples in the provided file to stdout.\n";
    out << " The tool also skips all information in the INFO column, but recalculates and sets allele count (AC) and allele number (AN) appropriately.\n";
    out << " Further, genotypes (GT) are expected to be the first entry in each sample column. Here, all information is kept." << endl;
    out << " Multi-allelics are not supported. All alleles differing from '0' are counted for AC." << endl;
    out << " Information will be printed to stderr." << endl;

    printVersion(out);
}

/* static */
void RemoveArgs::printVersion(ostream &out) {
    out << "This is version 0.1." << endl;
}
