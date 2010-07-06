// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#include "line_reader.hpp"
using namespace std;
using namespace Permory;

// Line_reader implementation
// ============================================================================
Line_reader::Line_reader(const char* fn, bool x) 
    : fn_(fn), extract_(x) 
{ 
    try {
        this->init();
    }
    catch (File_exception& e) {
        cerr << fn_ << ": " << e.what() << std::endl; 
        exit(-1);
    }
    catch (exception& e) {
        cerr << e.what() << std::endl;
        exit(-1);
    }
}

void Line_reader::init() throw (File_exception)
{
    buf_.resize(BUFFSIZE);
    if (!in_.empty()) {
        try {
            in_.reset(); //clear
        }
        catch (exception& e) {
            cerr << e.what(); exit(-1);
        }
    }
    // Check if file exists and extract file size
    File_info fi(fn_);
    if (fi.file_exists()) { 
        fileSize_ = fi.file_size(); 
    }
    else {
        //throw File_exception((string(fn_) + ": no such file." ).c_str());
        throw File_exception("no such file.");
    }

    // Create filter-device chain for the input stream
    if (extract_) in_.push(bio::zlib_decompressor()); //TODO
    //in_.push(Delim_input_filter(' '));
    in_.push(bio::file_source(fn_)); 
    assert (in_.is_complete());
    bool isOpen = in_.component<bio::file_source>(in_.size()-1)->is_open();
    if (!isOpen) 
        throw File_exception("unable to open file.");

    ++(*this); //init with first line
}

Line_reader& Line_reader::operator++() 
{
    nchar_ = 0;
    std::vector<char> v(0); //used in case line length exceeds buffer size

    while (!in_.eof()) {
        in_.getline(&*buf_.begin(), buf_.size());
        nchar_ += in_.gcount();

        if (in_.fail() && !in_.eof()) 
        {
            //FIXME undersize buffer gave segfault on MaRC (maybe buf_swap below?)
            // If we get here, line length exceeds buffer size: append copy
            // of what we have read so far to temporary buffer 'v ', and read on
            std::copy(buf_.begin(), buf_.begin() + strlen(&*buf_.begin()), 
                    std::back_inserter(v));
            in_.clear();
        }
        else
            break;
    }

    bool wasUndersized = v.size() > 0;
    if (wasUndersized) {
        // append rest to temporary buffer and then put everything back into 
        // "original" buffer, which also resizes the buffer for future use
        std::copy(buf_.begin(), buf_.begin() + strlen(&*buf_.begin()), 
                std::back_inserter(v));
        buf_.swap(v); 
    }
    ncharTotal_+= nchar_;
    lineCount_ += (nchar_ > 0);
    return *this;
}

Line_reader& Line_reader::operator+=(int n) 
{ 
    while(n-- > 0) 
        ++(*this); 
    return *this;
}

// ============================================================================
// Line_reader_alt implementation (not used)
// ============================================================================
Line_reader_alt::Line_reader_alt(const char* fn, bool x) 
    : fn_(fn), extract_(x) 
{ 
    try {
        init();
    }
    catch (File_exception& e) {
        cerr << e.what(); exit(-1);
    }
    catch (exception& e) {
        cerr << e.what(); exit(-1);
    }
}

void Line_reader_alt::init() throw (File_exception)
{
    buf_.resize(BUFFSIZE);
    if (!in_.empty()) {
        try {
            in_.reset(); //clear
        }
        catch (exception& e) {
            cerr << e.what(); exit(-1);
        }
    }
    // Check if file exists and extract file size
    File_info fi(fn_);
    if (fi.file_exists()) {
        fileSize_ = fi.file_size();
    }
    else {
        throw File_exception((string(fn_) + ": no such file." ).c_str());
    }

    // Create filter-device chain for the input stream
    //if (extract_) in_.push(bio::zlib_decompressor()); //TODO
    //in_.push(Delim_input_filter(' '));
    in_.push(bio::file_source(fn_)); 
    assert (in_.is_complete());
    bool isOpen = in_.component<bio::file_source>(in_.size()-1)->is_open();
    if (!isOpen) 
        throw File_exception((string(fn_) + ": unable to open file.").c_str());

    it_ = &in_;
    std::istreambuf_iterator<char> eos;
    eos_ = eos;
    ++(*this); //init with first line
}

Line_reader_alt& Line_reader_alt::operator++() 
{
    nchar_ = 0;
    std::vector<char>::iterator i = buf_.begin();
    while (it_ != eos_ && *it_ != '\n') {
        if (i == buf_.end()) {          //buffer full?
            //if buffer full:
            size_t sz = buf_.size();    //remember where we are
            buf_.resize(buf_.size() + BUFFSIZE); //enlarge buffer
            i = buf_.begin() + sz;      //read on from where we were
        }
        *i++ = *it_++;
        nchar_++;
    }
    *i = '\0'; //replace '\n'
    ncharTotal_+= nchar_;
    lineCount_ += (nchar_ > 0);

    if (it_ != eos_) {
        it_++; //point to next line
    }
    return *this;
}

Line_reader_alt& Line_reader_alt::operator+=(int n) 
{ 
    while(n-- > 0) 
        ++(*this); 
    return *this;
}

