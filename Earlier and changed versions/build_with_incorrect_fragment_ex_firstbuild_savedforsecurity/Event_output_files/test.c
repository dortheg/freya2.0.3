#include <experimental/string_view>
template<class _CharT, class _Traits = std::char_traits<_CharT> >                                                    
using basic_string_view = ::std::experimental::basic_string_view<_CharT,_Traits>;

int main() {                                                                     
   typedef basic_string_view<char> string_view;  
   return 0;
}