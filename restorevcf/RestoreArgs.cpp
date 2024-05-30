
#include <iostream>
#include <iomanip>

#include <boost/program_options.hpp>

#include "RestoreArgs.h"

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
/*static*/ RestoreArgs RestoreArgs::parseArgs(int argc, char *argv[]) {
    RestoreArgs args {argc, argv};

    if (args.count("help")) {
        args.printHelp(argv[0], cout);
        exit(EXIT_SUCCESS);
    }

    if (args.count("version")) {
        args.printVersion(cout);
        exit(EXIT_SUCCESS);
    }

    // set variables
    args.parseVars();

    return args;
}

ostream &operator<<(ostream &out, const RestoreArgs &args) {
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

RestoreArgs::RestoreArgs(int argc, char *argv[]) :
    opts_regular("Program options"),
	opts_hidden("Hidden options (only visible in debug mode)")
	{

    opts_regular.add_options()
    ("help,h", "produces this help message and exits")
    ("version,v", "prints version information and exits")
    ("fpass", "returns only variants with FILTER == PASS")
    ("rminfo", "removes all INFO fields (but still creates AF,AC,AN)")
    ("keepaa", "if present, the INFO field AAScore will be kept when removing the rest with --rminfo")
    ("macfilter", value<size_t>(&macfilter)->default_value(0), "only variants with a minor allele count >= value are returned")
    ("aafilter", value<float>(&aafilter)->default_value(0.0), "only variants with an AAScore >= value are returned")
    ("missfilter", value<float>(&missfilter)->default_value(0.0), "only variants with a genotype missingness rate < value are returned")
    ("filterunknown", "removes unknown alleles (named \"*\")")
    ("splitma", "splits multi-allelic variants into several bi-allelic ones, filling up with the reference '0'. Note, that this implies --rminfo.")
    ;

    opts_hidden.add_options()
    ("debug", "produce lots of debug output")
    ;

    parse(argc, argv);
}

void RestoreArgs::parse(int argc, char *argv[]) {
    bpo::options_description all_options;

    // combine all options
    all_options.add(opts_regular);
    all_options.add(opts_hidden);

    // do the actual parsing
    store(command_line_parser(argc, argv).options(all_options).run(), vars);
    notify(vars);

}

void RestoreArgs::parseVars() {

    if (vars.count("debug"))
        debug = true;

    // set bools
    if (vars.count("fpass"))
        fpass = true;
    if (vars.count("rminfo"))
        rminfo = true;
    if (vars.count("keepaa"))
        keepaa = true;
    if (vars.count("filterunknown"))
        filterunk = true;
    if (vars.count("splitma")) {
        splitma = true;
        rminfo = true; // this is automatically set when --splitma is used
    }
    // if we do not remove INFO, disable explicit keep of AAScore
    if (!rminfo)
        keepaa = false;

}

bool RestoreArgs::isDefined(const string &optname) const {
    bool found = false;
    found = !!this->opts_regular.find_nothrow(optname, false); // return null (-> false) if option has not been found
    found |= !!this->opts_hidden.find_nothrow(optname, false);
    return found;
}

void RestoreArgs::printHelp(const string &progname, ostream &out) const {
    out << "Usage: " << progname << " [options]" << endl << endl;
    out << opts_regular << endl;
    if (debug) {
        out << opts_hidden << endl;
    }
    out << endl;

    printVersion(out);
}

/* static */
void RestoreArgs::printVersion(ostream &out) {
    cout << "This is version 0.3." << endl;
}
