#include <string>
#include <cstdint>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <tuple>

#if defined(_MSC_VER) // needs to be first because msvc doesn't short-circuit after failing defined(__has_builtin)
#include <stdlib.h>
#  define bswap16(x)     _byteswap_ushort((x))
#  define bswap32(x)     _byteswap_ulong((x))
#  define bswap64(x)     _byteswap_uint64((x))
#elif (__GNUC__ > 4) || (__GNUC__ == 4 && __GNUC_MINOR__ >= 8)
#  define bswap16(x)     __builtin_bswap16((x))
#  define bswap32(x)     __builtin_bswap32((x))
#  define bswap64(x)     __builtin_bswap64((x))
#elif defined(__has_builtin) && __has_builtin(__builtin_bswap64)  /* for clang; gcc 5 fails on this and && shortcircuit fails; must be after GCC check */
#  define bswap16(x)     __builtin_bswap16((x))
#  define bswap32(x)     __builtin_bswap32((x))
#  define bswap64(x)     __builtin_bswap64((x))
#else
    /* even in this case, compilers often optimize by using native instructions */
    static inline uint16_t bswap16(uint16_t x) {
		return ((( x  >> 8 ) & 0xffu ) | (( x  & 0xffu ) << 8 ));
	}
    static inline uint32_t bswap32(uint32_t x) {
        return ((( x & 0xff000000u ) >> 24 ) |
                (( x & 0x00ff0000u ) >> 8  ) |
                (( x & 0x0000ff00u ) << 8  ) |
                (( x & 0x000000ffu ) << 24 ));
    }
    static inline uint64_t bswap64(uint64_t x) {
        return ((( x & 0xff00000000000000ull ) >> 56 ) |
                (( x & 0x00ff000000000000ull ) >> 40 ) |
                (( x & 0x0000ff0000000000ull ) >> 24 ) |
                (( x & 0x000000ff00000000ull ) >> 8  ) |
                (( x & 0x00000000ff000000ull ) << 8  ) |
                (( x & 0x0000000000ff0000ull ) << 24 ) |
                (( x & 0x000000000000ff00ull ) << 40 ) |
                (( x & 0x00000000000000ffull ) << 56 ));
    }
#endif

namespace py=pybind11;

std::tuple<py::array_t<float, py::array::c_style | py::array::forcecast>,
           py::array_t<float, py::array::c_style | py::array::forcecast>,
           py::array_t<float, py::array::c_style | py::array::forcecast>,
           py::array_t<float, py::array::c_style | py::array::forcecast>,
           py::array_t<float, py::array::c_style | py::array::forcecast>,
           py::array_t<float, py::array::c_style | py::array::forcecast>> unpack_rgb_sta (
    const char *raw_buffer,
    size_t sta_width,
    size_t sta_height,
    size_t sta_depth
) {

    // create all of the empty numpy arrays
    auto output_buffer_info = py::buffer_info(
            nullptr,            /* Pointer to data (nullptr -> ask NumPy to allocate!) */
            sizeof(float),     /* Size of one item */
            py::format_descriptor<float>::value, /* Buffer format */
            3,          /* How many dimensions? */
            {sta_depth, sta_width, sta_height},  /* Number of elements for each dimension */
            {sizeof(float) * sta_width * sta_height, sizeof(float) * sta_height, sizeof(float)}
            /* Strides for each dimension */
    );

    py::array_t <float> red_data = py::array_t<float>(output_buffer_info);
    py::buffer_info output_info = red_data.request();
    float *red_data_buffer = static_cast<float *> (output_info.ptr);

    py::array_t <float> red_error = py::array_t<float>(output_buffer_info);
    output_info = red_error.request();
    float *red_err_buffer = static_cast<float *> (output_info.ptr);

    py::array_t <float> green_data = py::array_t<float>(output_buffer_info);
    output_info = green_data.request();
    float *green_data_buffer = static_cast<float *> (output_info.ptr);

    py::array_t <float> green_error = py::array_t<float>(output_buffer_info);
    output_info = green_error.request();
    float *green_err_buffer = static_cast<float *> (output_info.ptr);

    py::array_t <float> blue_data = py::array_t<float>(output_buffer_info);
    output_info = blue_data.request();
    float *blue_data_buffer = static_cast<float *> (output_info.ptr);

    py::array_t <float> blue_error = py::array_t<float>(output_buffer_info);
    output_info = blue_error.request();
    float *blue_err_buffer = static_cast<float *> (output_info.ptr);

    size_t dwh_offset, wh_offset, write_offset;
    size_t read_offset = 3; // need to offset one 32 bit integer don't cares, and one 64 bit double don't care

    const uint32_t *raw_buffer_as_int = reinterpret_cast<const uint32_t *>(raw_buffer);
    uint32_t temp;
    for (size_t i = 0; i < sta_depth; ++i) {

	read_offset += 4; // need to offset two 32 bit integers and one 64 bit double that we don't care about

        dwh_offset = i * sta_width * sta_height;

        for (size_t j = 0; j < sta_width; ++j) {

            wh_offset = dwh_offset + j * sta_height;

            for (size_t k = 0; k < sta_height; ++k) {

                write_offset = wh_offset + k;

                temp = bswap32(*(raw_buffer_as_int+read_offset));
                *(red_data_buffer+write_offset) = *(reinterpret_cast<float *>(&temp));
                ++read_offset;

                temp = bswap32(*(raw_buffer_as_int+read_offset));
                *(red_err_buffer+write_offset) = *(reinterpret_cast<float *>(&temp));
                ++read_offset;

                temp = bswap32(*(raw_buffer_as_int+read_offset));
                *(green_data_buffer+write_offset) = *(reinterpret_cast<float *>(&temp));
                ++read_offset;

                temp = bswap32(*(raw_buffer_as_int+read_offset));
                *(green_err_buffer+write_offset) = *(reinterpret_cast<float *>(&temp));
                ++read_offset;

                temp = bswap32(*(raw_buffer_as_int+read_offset));
                *(blue_data_buffer+write_offset) = *(reinterpret_cast<float *>(&temp));
                ++read_offset;

                temp = bswap32(*(raw_buffer_as_int+read_offset));
                *(blue_err_buffer+write_offset) = *(reinterpret_cast<float *>(&temp));
                ++read_offset;
            }
        }
    }

    return std::make_tuple(red_data, red_error, green_data, green_error, blue_data, blue_error);
}
