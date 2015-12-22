package seq2c;

import java.util.concurrent.*;

public class Dispatcher {


    private static ExecutorService es;
    private static int threadsCount;

    public static void init(int threadsCount) {
        if (es != null) {
            throw new RuntimeException("Dispatcher is initialized");
        }

        es = new ThreadPoolExecutor(threadsCount, threadsCount,
                0L, TimeUnit.MILLISECONDS,
                new LinkedBlockingQueue<Runnable>());
        Dispatcher.threadsCount = threadsCount;
    }

    public static void shutdown() throws InterruptedException {
        if (es != null) {
            es.shutdown();
            es.awaitTermination(1, TimeUnit.DAYS);
            es = null;
        }
    }

    public static int getThreadsCount() {
        return threadsCount;
    }

    public static Service getService(int workers) {
        return new Service(workers);
    }

    public static class Service {


        private BlockingQueue<Worker> workers;

        private Service(int cnt) {
            this.workers = new LinkedBlockingQueue<Worker>(cnt);
        }


        public void submit(Runnable task) throws InterruptedException {
            Worker worker = new Worker(task);
            workers.put(worker);
            if (es == null) {
                worker.run();
            } else {
                es.submit(worker);
            }
        }

        public void await() throws InterruptedException, ExecutionException {
            for (Object worker : workers.toArray()) {
                ((Worker)worker).get();
            }
        }

        private class Worker extends FutureTask<Void> {

            Worker(Runnable task) {
                super(task, null);
            }

            @Override
            protected void done() {
                workers.remove(this);
            }
        }

    }

}
