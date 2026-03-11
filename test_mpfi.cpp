#include <mpfi.h>
#include <iostream>

int main() {
    // Test 1: Creating an mpfi via const mpfi_t&
    mpfi_t original;
    mpfi_init2(original, 64);
    mpfi_set_si(original, 5);
    
    // Test 2: Passing it to copy constructor
    mpfi_t copy;
    mpfi_init2(copy, mpfi_get_prec(original));
    mpfi_set(copy, original);
    
    // Test 3: What happens if we clear the local copy?
    mpfi_t temp;
    mpfi_init2(temp, mpfi_get_prec(copy));
    mpfi_add(temp, copy, original);
    // Now if we create an FMpfi from temp and then clear temp...
    // The pointer inside FMpfi.val would be pointing to freed memory
    mpfi_clear(temp);  // This would invalidate any reference to temp's data
    
    mpfi_clear(copy);
    mpfi_clear(original);
    return 0;
}
