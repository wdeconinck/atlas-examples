#include <filesystem>
#include <fstream>
#include <iostream>

#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/library.h"
#include "atlas/mesh.h"
#include "atlas/output/Gmsh.h"
#include "atlas/runtime/Log.h"

#include <stdio.h>
#include "eccodes.h"

#include "atlas/io/atlas-io.h"

// --------------------------------------------------------------------------------------------------------------

class GribFileReader {
public:
    GribFileReader(const std::string& path) {
        if (!std::filesystem::exists(path)) {
            std::cerr << path << " does not exist" << std::endl;
            std::exit(1);
        }
        file_ = fopen(path.c_str(), "rb");
        if (!file_) {
            std::cerr << "Could not open file " << path << std::endl;
            std::exit(1);
        }
        codes_count_in_file(nullptr, file_, &count_);
        int err;
        grib_ = codes_handle_new_from_file(nullptr, file_, PRODUCT_GRIB, &err);
        if (!grib_) {
            std::cerr << "Could not create grib handle. Error: " << err << std::endl;
            std::exit(1);
        }
    }
    ~GribFileReader() {
        if (grib_)
            codes_handle_delete(grib_);
        if (file_)
            fclose(file_);
    }
    size_t get_values_size() {
        size_t size;
        codes_get_size(grib_, "values", &size);
        return size;
    }
    void get_values(double values[], size_t size) {
        if (size < get_values_size()) {
            std::cerr << "values[] has been allocated too small." << std::endl;
            std::exit(1);
        }
        codes_get_double_array(grib_, "values", values, &size);
    }
    std::string get_gridname() {
        std::string gridType = get_string("gridType");
        if (gridType == "reduced_gg") {
            if (get_long("isOctahedral")) {
                return "O" + std::to_string(get_long("N"));
            }
            else {
                return "N" + std::to_string(get_long("N"));
            }
        }
        else if (gridType == "regular_gg") {
            return "F" + std::to_string(get_long("N"));
        }
        return gridType;
    }
    std::string get_string(const std::string& key) {
        char buffer[64];
        size_t size = sizeof(buffer);
        codes_get_string(grib_, key.c_str(), buffer, &size);
        return std::string(buffer);
    };
    long get_long(const std::string& key) {
        long value;
        codes_get_long(grib_, key.c_str(), &value);
        return value;
    };
    bool next_message() {
        if (index_ == count_) {
            return false;
        }
        if (grib_) {
            codes_handle_delete(grib_);
        }
        int err;
        grib_ = codes_handle_new_from_file(nullptr, file_, PRODUCT_GRIB, &err);
        if (!grib_) {
            std::cerr << "Could not create grib handle in next_message(). Error: " << err << std::endl;
            std::exit(1);
        }
        index_++;
        return true;
    }
    int count() const { return count_; }

private:
    FILE* file_         = nullptr;
    codes_handle* grib_ = nullptr;

    int count_{0};
    int index_{1};
};

// --------------------------------------------------------------------------------------------------------------

void convert_grib_to_atlas_io(const std::string& grib_file, const std::string& atlas_io_file) {
    if (atlas::mpi::rank() == 0) {
        GribFileReader grib{grib_file};
        size_t nfld   = grib.count();
        auto gridname = grib.get_gridname();
        atlas::io::RecordWriter atlas_io_writer;
        auto compression = atlas::util::Config("compression", "lz4");
        atlas_io_writer.set("grid.name", atlas::io::ref(gridname), compression);
        std::cout << "grid.name: " << gridname << std::endl;
        std::cout << "grid.size: " << grib.get_values_size() << std::endl;

        atlas_io_writer.set("fields.size", nfld);
        std::cout << "fields: " << std::endl;

        std::vector<double> values(grib.get_values_size());

        for (size_t jfld = 0; jfld < nfld; ++jfld) {
            grib.get_values(values.data(), values.size());
            auto fieldname = grib.get_string("shortName");
            auto fielddesc = grib.get_string("name");
            long level     = grib.get_long("level");

            std::cout << "    " << std::setw(5) << std::left << jfld << std::setw(16) << std::left << fieldname << " ["
                      << level << "] " << fielddesc << std::endl;

            std::string field = "fields[" + std::to_string(jfld) + "]";
            atlas_io_writer.set(field + ".name", fieldname);
            atlas_io_writer.set(field + ".description", fielddesc);
            atlas_io_writer.set(field + ".level", level);
            atlas_io_writer.set(field + ".array", atlas::io::copy(values), compression);
            grib.next_message();
        }

        atlas_io_writer.write(atlas_io_file);
    }
    atlas::mpi::comm().barrier();
}

// --------------------------------------------------------------------------------------------------------------

class AtlasIOFileReader {
public:
    AtlasIOFileReader(const std::string& path): record_(path) {}

    template <typename T>
    T read(const std::string& key) {
        T value;
        record_.read(key, value).wait();
        return value;
    }

    template <typename T>
    void read(const std::string& key, T& value) {
        record_.read(key, value).wait();
    }

    atlas::io::RecordReader& record() { return record_; }

private:
    atlas::io::RecordReader record_;
};

// --------------------------------------------------------------------------------------------------------------

struct CommandLineOptions {
    bool gmsh                    = false;
    std::string gmsh_file        = "out.msh";
    std::string gmsh_coordinates = "xy";
    std::string atlas_io_file    = "out.atlas";
    std::string grib_file        = "in.grib";
    CommandLineOptions(int argc, char* argv[]) {
        int c;
        auto get_opt = [&](const std::string_view& opt, std::string& value) {
            std::string_view arg(argv[c]);
            if (arg == opt) {
                if (c < argc - 1) {
                    if (arg[0] != '-') {
                        value = argv[++c];
                    }
                }
            }
        };
        for (c = 1; c < argc; ++c) {
            if (c == 1) {
                grib_file = argv[c];
                continue;
            }
            if (std::string_view(argv[c]) == "--gmsh") {
                gmsh = true;
            }
            get_opt("-o", atlas_io_file);
            get_opt("--output", atlas_io_file);
            get_opt("--gmsh", gmsh_file);
            get_opt("--coordinates", gmsh_coordinates);
        }
    }
};

// --------------------------------------------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    atlas::initialize(argc, argv);

    CommandLineOptions args(argc, argv);

    convert_grib_to_atlas_io(args.grib_file, args.atlas_io_file);

    if (args.gmsh) {
        atlas::Log::info() << "Output to gmsh file " << args.gmsh_file << std::endl;
        // For verification, and for example use, we can load the fields
        // from the generated file and visualize them with gmsh

        AtlasIOFileReader reader{args.atlas_io_file};

        auto grid = atlas::Grid(reader.read<std::string>("grid.name"));
        auto mesh = atlas::Mesh(grid);

        auto gmsh = atlas::output::Gmsh(args.gmsh_file, atlas::util::Config("coordinates", args.gmsh_coordinates));
        gmsh.write(mesh);

        auto fs = atlas::functionspace::NodeColumns(mesh);

        atlas::FieldSet fields;

        size_t field_size = reader.read<size_t>("fields.size");
        for (size_t i = 0; i < field_size; ++i) {
            std::string prefix = "fields[" + std::to_string(i) + "]";
            long level         = reader.read<long>(prefix + ".level");
            std::string name   = reader.read<std::string>(prefix + ".name") + "[" + std::to_string(level) + "]";
            auto field         = fs.createField<double>(atlas::option::name(name));
            {
                auto field_global = fs.createField(field, atlas::option::global());
                if (atlas::mpi::rank() == 0) {
                    reader.read(prefix + ".array", field_global.array());
                }
                fs.scatter(field_global, field);
            }
            fields.add(field);
        }

        fs.haloExchange(fields);
        gmsh.write(fields);
    }

    atlas::finalize();
    return 0;
}

// --------------------------------------------------------------------------------------------------------------
