# Contributing Guidelines

Thank you for your interest in contributing to our project. Whether it's a bug report, new feature, correction, or additional
documentation, we greatly value feedback and contributions from our community.

Please read through this document before submitting any issues or pull requests to ensure we have all the necessary
information to effectively respond to your bug report or contribution.


## Adding additional modules

Follow these steps to add additional modules

1. Create a new folder under `infrastructure/docker` with your Dockerfile and any supporting files. Look at the `run.sh` helper script included with other modules for and example of how to pass data back and forth with Amazon S3.
2. Copy one of the `batch-protein-folding-cdn-module-*.yaml` files in `/infrastructure/cloudformation` and name it after your new module. Update the resource names, CodeBuild details, and Job Definition as needed.
3. If your modeul needs data on the FSx for Lustre file system (e.g. model weights or reference data), create a new download script in `/infrastructure/docker/download/script`. You'll also need to add the name of the script to the download function defined in `prep_databases.py` and `infrastructure/cloudformation/batch-protein-folding-cfn-download.yaml`.
4. Create a new module for submitting your job in `src/batchfold` and tests in `tests`.
5. Build and deploy your stack using the following AWS CLI commands:
```
aws cloudformation package --template-file infrastructure/cloudformation/batch-protein-folding-cfn-root.yaml --region <YOUR REGION> --output-template infrastructure/cloudformation/batch-protein-folding-cfn-packaged.yaml --s3-bucket <YOUR S3 BUCKET NAME>

aws cloudformation deploy --template-file infrastructure/cloudformation/batch-protein-folding-cfn-packaged.yaml --region <YOUR REGION> --capabilities CAPABILITY_IAM --stack-name <YOUR STACK NAME>
```
6. Find the URL for your Code Commit repository, add it as a remote to your local git repository, and push your local changes. This will overwrite the "public" code with your updates.
7. Start the CodeBuild project associated to your new module to create the container.

## Reporting Bugs/Feature Requests

We welcome you to use the GitHub issue tracker to report bugs or suggest features.

When filing an issue, please check existing open, or recently closed, issues to make sure somebody else hasn't already
reported the issue. Please try to include as much information as you can. Details like these are incredibly useful:

* A reproducible test case or series of steps
* The version of our code being used
* Any modifications you've made relevant to the bug
* Anything unusual about your environment or deployment


## Contributing via Pull Requests
Contributions via pull requests are much appreciated. Before sending us a pull request, please ensure that:

1. You are working against the latest source on the *main* branch.
2. You check existing open, and recently merged, pull requests to make sure someone else hasn't addressed the problem already.
3. You open an issue to discuss any significant work - we would hate for your time to be wasted.

To send us a pull request, please:

1. Fork the repository.
2. Modify the source; please focus on the specific change you are contributing. If you also reformat all the code, it will be hard for us to focus on your change.
3. Ensure local tests pass.
4. Commit to your fork using clear commit messages.
5. Send us a pull request, answering any default questions in the pull request interface.
6. Pay attention to any automated CI failures reported in the pull request, and stay involved in the conversation.

GitHub provides additional document on [forking a repository](https://help.github.com/articles/fork-a-repo/) and
[creating a pull request](https://help.github.com/articles/creating-a-pull-request/).


## Finding contributions to work on
Looking at the existing issues is a great way to find something to contribute on. As our projects, by default, use the default GitHub issue labels (enhancement/bug/duplicate/help wanted/invalid/question/wontfix), looking at any 'help wanted' issues is a great place to start.


## Code of Conduct
This project has adopted the [Amazon Open Source Code of Conduct](https://aws.github.io/code-of-conduct).
For more information see the [Code of Conduct FAQ](https://aws.github.io/code-of-conduct-faq) or contact
opensource-codeofconduct@amazon.com with any additional questions or comments.


## Security issue notifications
If you discover a potential security issue in this project we ask that you notify AWS/Amazon Security via our [vulnerability reporting page](http://aws.amazon.com/security/vulnerability-reporting/). Please do **not** create a public github issue.


## Licensing

See the [LICENSE](LICENSE) file for our project's licensing. We will ask you to confirm the licensing of your contribution.
