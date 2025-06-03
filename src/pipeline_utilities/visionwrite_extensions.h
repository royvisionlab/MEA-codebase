#include <string>
#include <cstdint>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>


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

//Start of mike adding
// #ifdef _MSC_VER
//
// #include <stdlib.h>
// // #define bswap_32(x) _byteswap_ulong(x)
// // #define bswap_64(x) _byteswap_uint64(x)
// #define __builtin_bswap_32(x) _byteswap_ulong(x)
// #define __builtin_bswap_64(x) _byteswap_uint64(x)
//
//
// #elif defined(__APPLE__)
//
// // Mac OS X / Darwin features
// #include <libkern/OSByteOrder.h>
// #define bswap_32(x) OSSwapInt32(x)
// #define bswap_64(x) OSSwapInt64(x)
//
// #elif defined(__sun) || defined(sun)
//
// #include <sys/byteorder.h>
// #define bswap_32(x) BSWAP_32(x)
// #define bswap_64(x) BSWAP_64(x)
//
// #elif defined(__FreeBSD__)
//
// #include <sys/endian.h>
// #define bswap_32(x) bswap32(x)
// #define bswap_64(x) bswap64(x)
//
// #elif defined(__OpenBSD__)
//
// #include <sys/types.h>
// #define bswap_32(x) swap32(x)
// #define bswap_64(x) swap64(x)
//
// #elif defined(__NetBSD__)
//
// #include <sys/types.h>
// #include <machine/bswap.h>
// #if defined(__BSWAP_RENAME) && !defined(__bswap_32)
// #define bswap_32(x) bswap32(x)
// #define bswap_64(x) bswap64(x)
// #endif
//
// #else
//
// #include <byteswap.h>
//
// #endif
// End of Mike Adding

py::bytes pack_sta_buffer_color (
    py::array_t<float, py::array::c_style | py::array::forcecast>& red_sta,
    py::array_t<float, py::array::c_style | py::array::forcecast>& red_err,
    py::array_t<float, py::array::c_style | py::array::forcecast>& green_sta,
    py::array_t<float, py::array::c_style | py::array::forcecast>& green_err,
    py::array_t<float, py::array::c_style | py::array::forcecast>& blue_sta,
    py::array_t<float, py::array::c_style | py::array::forcecast>& blue_err,
    double refresh_time) {

    /*
    Note that we require that the arrrays have a different array order than the STAs
        returned by visionloader. This is because the original order is a nightmare
        for cache in memory accesses (and hence comically slow for fine
        stimuli)

    We will require an array order swap before using this function
    */

    py::buffer_info red_sta_info = red_sta.request();
    float *red_data_ptr = static_cast<float *> (red_sta_info.ptr);

    size_t sta_depth = red_sta_info.shape[0];
    size_t sta_width = red_sta_info.shape[1];
    size_t sta_height = red_sta_info.shape[2];

    size_t n_output_entries = 6 * sta_width * sta_height * sta_depth + sta_depth * 4;


    py::buffer_info red_err_info = red_err.request();
    float *red_err_ptr = static_cast<float *> (red_err_info.ptr);

    py::buffer_info green_sta_info = green_sta.request();
    float *green_data_ptr = static_cast<float *> (green_sta_info.ptr);

    py::buffer_info green_err_info = green_err.request();
    float *green_err_ptr = static_cast<float *> (green_err_info.ptr);

    py::buffer_info blue_sta_info = blue_sta.request();
    float *blue_data_ptr = static_cast<float *> (blue_sta_info.ptr);

    py::buffer_info blue_err_info = blue_err.request();
    float *blue_err_ptr = static_cast<float *> (blue_err_info.ptr);

    uint32_t *output_buffer = new uint32_t[n_output_entries];

    size_t depth_offset, width_depth_offset, read_offset;
    uint64_t refresh_temp;
    size_t write_idx = 0;
    for (size_t i = 0; i < sta_depth; ++i) {

        depth_offset = i * (sta_width * sta_height);

        output_buffer[write_idx++] = bswap32(static_cast<uint32_t>(sta_width));
        output_buffer[write_idx++] = bswap32(static_cast<uint32_t>(sta_height));

        refresh_temp = bswap64(*(reinterpret_cast<uint64_t *>(&refresh_time)));
        output_buffer[write_idx++] = static_cast<uint32_t> (refresh_temp >> 32);
        output_buffer[write_idx++] = static_cast<uint32_t> (refresh_temp & 0xFFFF);

        for (size_t j = 0; j < sta_width; ++j)  {

            width_depth_offset = j * sta_height + depth_offset;

            for (size_t k = 0; k < sta_height; ++k) {

                read_offset = width_depth_offset + k;

                output_buffer[write_idx++] = bswap32(*(reinterpret_cast<uint32_t *> (red_data_ptr + read_offset)));
                output_buffer[write_idx++] = bswap32(*(reinterpret_cast<uint32_t *> (red_err_ptr + read_offset)));

                output_buffer[write_idx++] = bswap32(*(reinterpret_cast<uint32_t *> (green_data_ptr + read_offset)));
                output_buffer[write_idx++] = bswap32(*(reinterpret_cast<uint32_t *> (green_err_ptr + read_offset)));

                output_buffer[write_idx++] = bswap32(*(reinterpret_cast<uint32_t *> (blue_data_ptr + read_offset)));
                output_buffer[write_idx++] = bswap32(*(reinterpret_cast<uint32_t *> (blue_err_ptr + read_offset)));

            }
        }
    }


    return py::bytes(reinterpret_cast<char *> (output_buffer), n_output_entries * 4);
}
