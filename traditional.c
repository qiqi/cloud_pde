/* A simple server in the internet domain using TCP
   The port number is passed as an argument */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include <sys/types.h> 
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h> 
#include <pthread.h>

void error(const char *msg)
{
    perror(msg);
    exit(1);
}

struct start_server_args {
    struct hostent * server_L;
    struct hostent * server_R;
    int sockfd_L;
    int sockfd_R;
};

void start_server(struct start_server_args * args)
{
    int sockfd = socket(AF_INET, SOCK_STREAM, 0);
    assert(sockfd >= 0);

    struct sockaddr_in serv_addr;
    bzero((char *) &serv_addr, sizeof(serv_addr));
    serv_addr.sin_family = AF_INET;
    serv_addr.sin_addr.s_addr = INADDR_ANY;
    serv_addr.sin_port = htons(30000);
    assert(bind(sockfd, (struct sockaddr *) &serv_addr, sizeof(serv_addr))
           >= 0);
    listen(sockfd,10);

    for (int i_sock = 0; i_sock < 2; ++ i_sock) {
        struct sockaddr_in cli_addr;
        socklen_t clilen = sizeof(cli_addr);
        int newsockfd = accept(sockfd, (struct sockaddr*)&cli_addr, &clilen);
        assert(newsockfd >= 0);
        if (args->server_L && memcmp(&cli_addr.sin_addr.s_addr,
                                     args->server_L->h_addr_list[0],
                                     args->server_L->h_length) == 0) {
            args->sockfd_L = newsockfd;
        } else {
            assert(args->server_R && memcmp(&cli_addr.sin_addr.s_addr,
                                            args->server_R->h_addr_list[0],
                                            args->server_R->h_length) == 0);
            args->sockfd_R = newsockfd;
        }
    }
    pthread_exit(0);
}

int connect_to(struct hostent * server)
{
    int sockfd = -1;
    if (server) {
        struct sockaddr_in serv_addr;
        bzero((char *) &serv_addr, sizeof(serv_addr));
        serv_addr.sin_family = AF_INET;
        bcopy(server->h_addr_list[0], &serv_addr.sin_addr.s_addr,
              server->h_length);
        serv_addr.sin_port = htons(8000);
        while (sockfd < 0) {
            connect(sockfd, (struct sockaddr*)&serv_addr,
                              sizeof(serv_addr));
        }
    }
    return sockfd;
}

int main(int argc, char *argv[])
{
    assert (argc == 3);
    struct hostent * server_L = gethostbyname(argv[1]);
    struct hostent * server_R = gethostbyname(argv[2]);

    struct start_server_args receiver;
    receiver.server_L = server_L;
    receiver.server_R = server_R;
    receiver.sockfd_L = receiver.sockfd_R = -1;
    pthread_t server_thread;
    pthread_create(&server_thread, 0, (void(*))&start_server, &receiver);

    int sockfd_L = connect_to(server_L);
    int sockfd_R = connect_to(server_R);

    while(1) {
        if (server_L && receiver.sockfd_L <0) {}
        if (server_R && receiver.sockfd_R <0) {}
        else { break; }
    }
    printf("L: (%d %d); R: (%d %d)", receiver.sockfd_L, sockfd_L,
                                     receiver.sockfd_R, sockfd_R);

    return 0; 
}
