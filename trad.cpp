#include<iostream>
#include<zmq.hpp>

int main()
{
    zmq::context_t ctx(1);
    zmq::socket_t socket_L(ctx, ZMQ_REQ);
    zmq::socket_t socket_R(ctx, ZMQ_REP);
    socket_L.connect("tcp://128.0.0.1:5555");
    socket_R.bind("tcp://eth0:5555");
    zmq::message_t msg(100);
    memset(msg.data(), 0, 100);
    socket_L.send(msg);
    return 0;
}
