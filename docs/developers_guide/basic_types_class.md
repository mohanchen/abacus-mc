# Basic Tool Classes Design Guide for ABACUS

## Overview

This document provides guidelines for designing and implementing basic tool classes in the ABACUS codebase, focusing on best practices for memory management, code style, and testing. These guidelines apply to all basic mathematical and utility classes, including but not limited to:

- vector3.h
- matrix.h
- timer.h
- ndarray.h
- realarray.h
- complexarray.h
- complexmatrix.h
- matrix3.h
- intarray.h
- formatter.h
- math_chebyshev.h

While this guide uses `IntArray` as an example for illustration purposes, the principles and practices described here are applicable to all basic tool classes in ABACUS.

## Memory Management

### 1. Exception Handling for Memory Allocation

Always use try-catch blocks when allocating memory to handle `std::bad_alloc` exceptions gracefully:

```cpp
try
{
    ptr = new int[size];
    zero_out();
}
catch (const std::bad_alloc& e)
{
    std::cerr << "Allocation error for IntArray: " << e.what() << std::endl;
    ptr = nullptr;
    size = 0;
    throw;
}
assert(ptr != nullptr);
```

### 2. Two-Stage Memory Allocation

When reallocating memory (e.g., in `create` methods), use a two-stage approach to ensure that the original object remains valid if memory allocation fails:

```cpp
int* new_ptr = nullptr;
try
{
    new_ptr = new int[size];
}
catch (const std::bad_alloc& e)
{
    std::cerr << "Allocation error in IntArray::create: " << e.what() << std::endl;
    assert(new_ptr != nullptr);
    return;
}
delete[] ptr;
ptr = new_ptr;
zero_out();
```

### 3. Null Pointer Checks

Always check for null pointers before accessing memory, especially in methods that might be called on objects with failed memory allocation:

```cpp
void IntArray::zero_out()
{
    if (size <= 0 || ptr == nullptr)
    {
        return;
    }
    for (int i = 0; i < size; i++)
    {
        ptr[i] = 0;
    }
    return;
}
```

## Class Design

### 1. Copy Constructor

Implement a copy constructor to avoid shallow copy issues:

```cpp
IntArray::IntArray(const IntArray& other)
{
    size = other.size;
    dim = other.dim;
    bound1 = other.bound1;
    bound2 = other.bound2;
    bound3 = other.bound3;
    bound4 = other.bound4;
    bound5 = other.bound5;
    bound6 = other.bound6;
    try
    {
        ptr = new int[size];
        for (int i = 0; i < size; i++)
        {
            ptr[i] = other.ptr[i];
        }
    }
    catch (const std::bad_alloc& e)
    {
        std::cerr << "Allocation error in IntArray copy constructor: " << e.what() << std::endl;
        ptr = nullptr;
        size = 0;
        throw;
    }
    assert(ptr != nullptr);
}
```

### 2. Move Semantics

Implement move constructor and move assignment operator to improve performance:

```cpp
// Move constructor
IntArray::IntArray(IntArray&& other) noexcept
    : size(other.size),
      dim(other.dim),
      bound1(other.bound1),
      bound2(other.bound2),
      bound3(other.bound3),
      bound4(other.bound4),
      bound5(other.bound5),
      bound6(other.bound6),
      ptr(other.ptr)
{
    other.ptr = nullptr;
    other.size = 0;
    other.dim = 0;
    other.bound1 = other.bound2 = other.bound3 = other.bound4 = other.bound5 = other.bound6 = 0;
}

// Move assignment operator
IntArray& IntArray::operator=(IntArray&& other) noexcept
{
    if (this != &other)
    {
        freemem();
        size = other.size;
        dim = other.dim;
        bound1 = other.bound1;
        bound2 = other.bound2;
        bound3 = other.bound3;
        bound4 = other.bound4;
        bound5 = other.bound5;
        bound6 = other.bound6;
        ptr = other.ptr;
        other.ptr = nullptr;
        other.size = 0;
        other.dim = 0;
        other.bound1 = other.bound2 = other.bound3 = other.bound4 = other.bound5 = other.bound6 = 0;
    }
    return *this;
}
```

### 3. Boundary Checks

Add boundary checks to prevent out-of-bounds access:

```cpp
int& IntArray::operator()(const int d1, const int d2)
{
    assert(d1 >= 0 && d1 < bound1);
    assert(d2 >= 0 && d2 < bound2);
    return ptr[d1 * bound2 + d2];
}
```

## Code Style

### 1. Brace Style

Use separate lines for braces:

```cpp
if (condition)
{
    // code here
}

void function()
{
    // code here
}
```

### 2. Indentation

Use spaces instead of tabs for indentation (4 spaces per indent level).

### 3. Comments

Use English for comments and document important functionality. Follow Doxygen-style documentation for classes and methods.

## Testing

### 1. Unit Tests

Write comprehensive unit tests for all classes, including:
- Constructor tests
- Method tests
- Exception handling tests
- Edge case tests

### 2. Test Class Initialization

Use constructor initialization lists for test classes to improve compatibility:

```cpp
class IntArrayTest : public testing::Test
{
protected:
    ModuleBase::IntArray a2, a3, a4, a5, a6;
    int aa;
    int bb;
    int count0;
    int count1;
    const int zero;

    IntArrayTest() : aa(11), bb(1), zero(0)
    {
    }
};
```

## Best Practices

1. **Single Responsibility Principle**: Each class should have a single, well-defined responsibility.
2. **Encapsulation**: Hide implementation details and expose only necessary interfaces.
3. **Error Handling**: Handle errors gracefully, especially memory allocation failures.
4. **Performance**: Use move semantics and other performance optimizations where appropriate.
5. **Testing**: Write comprehensive tests for all functionality.
6. **Code Style**: Follow consistent code style guidelines.
7. **Documentation**: Document classes and methods to improve maintainability.
8. **Compatibility**: Ensure code is compatible with C++11 standard.
9. **Portability**: Write code that works across different platforms.
10. **Reusability**: Design classes to be reusable in different contexts.

## Application to Other Basic Tool Classes

While this guide uses `IntArray` as an example, these principles apply to all basic tool classes in ABACUS. For example:

- **vector3.h**: Apply the same memory management and error handling principles, with additional focus on vector operations and operator overloading.
- **matrix.h**: Extend the memory management practices to 2D arrays, with additional considerations for matrix operations.
- **timer.h**: Focus on static member management and time measurement accuracy.
- **ndarray.h**: Apply the same principles to multi-dimensional arrays, with additional considerations for shape manipulation.
- **formatter.h**: Focus on string manipulation and formatting, with attention to performance and usability.
- **math_chebyshev.h**: Apply the principles to template classes, with additional focus on mathematical algorithm implementation.

By following these guidelines, you can ensure that all basic tool classes in ABACUS are well-designed, robust, and maintainable.
