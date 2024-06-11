
# Package infrastructure for deployment
aws cloudformation package --template-file infrastructure/cloudformation/batch-protein-folding-cfn-root.yaml --s3-bucket archon-iac --s3-prefix main --output-template-file ../protein-folding.cfn.yaml --profile archon

# Update stack
aws cloudformation update-stack --stack-name archon-iac-galaxy --template-body file://../protein-folding.cfn.yaml --capabilities CAPABILITY_NAMED_IAM --profile archon --disable-rollback --output json --profile archon