
#ifndef ARGS_H_
#define ARGS_H_

#include <string>
#include <vector>

#include <boost/program_options.hpp>

namespace bpo = boost::program_options;

using namespace std;

/**
 * Class for storing and retrieving command-line arguments.
 */
class RestoreArgs {
public:

	static RestoreArgs parseArgs(int argc, char *argv[]);

    /**
     * Returns the argument with the given (long) name. The template argument
     * specifies the return type you want the parameter casted to.
     * @param name argument name
     * @return the variable value
     */
    template<typename T> T get(const std::string &name) const {
        auto where = vars.find(name);
        if(where == std::end(vars)) {
            if(!isDefined(name))
                throw std::invalid_argument("Option undefined: " + name + " (This is a bug)");
            else
                throw std::out_of_range("Option has not been specified and does not have a default value associated: " + name);
        }
        return where->second.as<T>();
    }

    /**
     * Counts the number of argument occurences. Mainly useful for boolean switches
     * and/or counting flags, such as verbosity or debug mode.
     * @param name argument name
     * @return argument value
     */
    unsigned int count(const std::string &name) const {
        if(!isDefined(name))
            throw std::invalid_argument("Option undefined: " + name + " (This is a bug)");
        return vars.count(name);
    }

    bool operator()(const std::string &name) const {
        return count(name) > 0;
    }

    /**
     * Prints a help message.
     * @param progname the program name, usually argv[0]
     * @param out output stream, e.g. cout
     */
    void printHelp(const std::string &progname, std::ostream &out) const;
    static void printVersion(std::ostream &out);

    RestoreArgs(RestoreArgs&& other) = default;

    bool fpass = false;
    bool rminfo = false;
    bool keepaa = false;
    size_t macfilter = 0;
    float aafilter = 0;
    float missfilter = 0;
    bool filterunk = false;
    bool splitma = false;

    bool debug = false;

protected:
    /** Constructs the arguments list and adds all defined options */
    RestoreArgs();
    RestoreArgs(RestoreArgs const &);
    void operator=(RestoreArgs const &);

    void parse(int argc, char *argv[]);
    void parseVars();
    bool isDefined(const std::string &optname) const;

    bpo::options_description opts_regular;        /**< regular options, shown on help output */
    bpo::options_description opts_hidden;         /**< hidden options */

    bpo::variables_map vars;    /**< this is where the option values go */

    /** dump all options to the given stream */
    friend std::ostream &operator<<(std::ostream &out, const RestoreArgs &args);

private:
    /** parses the main() options into the variable map */
    RestoreArgs(int argc, char *argv[]);

};

#endif /* ARGS_H_ */
