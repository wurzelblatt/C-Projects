#include <iostream>

using namespace std;

class Number {
	private:
	int x;
	public:
	Number(int x) {
		this->x = x;
	}
	void* operator new (size_t size) {
			void *ptr = ::new int[size]; //Using global new operator
			cout<<"Memory allocated of size: "<<size<<endl;
			return ptr;
	}
	void operator delete(void *ptr) {
		::delete(ptr);   //Using global delete operator
		cout<<"Memory deallocated"<<endl;
	}
	void display() {
		cout<<"x = "<<x<<endl;
	}
};
int main()
{
	Number *n = new  Number(10);  //Invokes overloaded new operator
	n->display();
	delete n; 	//Invokes overloaded delete operator
	return 0;
}
